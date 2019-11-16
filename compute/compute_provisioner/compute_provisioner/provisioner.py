import datetime
import json
import logging
import threading
import time
from collections import OrderedDict
from uuid import uuid4

import redis
import yaml
import zmq
from kubernetes import client
from kubernetes import config

from compute_provisioner import constants
from compute_provisioner import kube_cluster

logger = logging.getLogger('provisioner')

CLEANUP_PHASES = ('Succeeded', 'Failed')


class ResourceTimeoutError(Exception):
    def __init__(self, name):
        super(ResourceTimeoutError, self).__init__('Timed out waitin for resources {!s}'.format(name))


class ResourceAllocationError(Exception):
    def __init__(self, name):
        self.name = name

        super(ResourceAllocationError, self).__init__('Error allocating resource {!s}'.format(name))


def json_encoder(x):
    def default(y):
        return {'data': y.decode(constants.ENCODING)}

    return json.dumps(x, default=default)


def json_decoder(x):
    def object_hook(y):
        if 'data' in y:
            return y['data'].encode(constants.ENCODING)

        return y

    return json.loads(x, object_hook=object_hook)


class WaitingAck(object):
    def __init__(self, version, frames):
        self.version = version
        self.frames = frames
        self.expiry = time.time() + constants.ACK_TIMEOUT


class Worker(object):
    def __init__(self, address, version):
        self.address = address
        self.version = version
        self.expiry = time.time() + constants.HEARTBEAT_INTERVAL * constants.HEARTBEAT_LIVENESS


class WorkerQueue(object):
    def __init__(self):
        self.queue = OrderedDict()

        self.heartbeat_at = time.time() + constants.HEARTBEAT_INTERVAL

    def ready(self, worker):
        self.queue.pop(worker.address, None)

        self.queue[worker.address] = worker

        logger.debug('Setting worker %r to expire at %r', worker.address, worker.expiry)

    def purge(self):
        t = time.time()

        expired = []

        for address, worker in self.queue.items():
            if t > worker.expiry:
                expired.append(address)

        for address in expired:
            self.queue.pop(address, None)

            logger.info('Idle worker %r expired', address)

    @property
    def heartbeat_time(self):
        return time.time() >= self.heartbeat_at

    def heartbeat(self, socket):
        for address in self.queue:
            socket.send_multipart([address, constants.HEARTBEAT])

            logger.debug('Sent heartbeat to %r', address)

        self.heartbeat_at = time.time() + constants.HEARTBEAT_INTERVAL

        logger.debug('Setting next heartbeat at %r', self.heartbeat_at)

    def next(self, version):
        candidates = [x.address for x in self.queue.values() if x.version == version]

        logger.debug('Candidates for next worker handling version %r: %r', version, candidates)

        address = candidates.pop(0)

        self.queue.pop(address, None)

        return address


class Provisioner(threading.Thread):
    def __init__(self, frontend_port, backend_port, redis_host, namespace, lifetime, resource_timeout, wait_for_resources, **kwargs):
        super(Provisioner, self).__init__(target=self.monitor, args=(frontend_port, backend_port, redis_host))
        self.wait_for_resources = wait_for_resources

        self.namespace = namespace

        self.lifetime = lifetime

        self.resource_timeout = resource_timeout

        self.context = None

        self.frontend = None

        self.backend = None

        self.poll_both = None

        self.workers = WorkerQueue()

        self.redis = None

        self.waiting_ack = {}

        config.load_incluster_config()

        self.core = client.CoreV1Api()

        self.apps = client.AppsV1Api()

        self.extensions = client.ExtensionsV1beta1Api()

    def initialize(self, frontend_port, backend_port, redis_host):
        """ Initializes the load balancer.

        We initialize the zmq context, create the frontend/backend sockets and
        creates their respective pollers.

        Args:
            frontend_port: An int port number used to listen for frontend
                connections.
            backend_port: An int port number used to listen for backend
                connections.
        """
        self.context = zmq.Context(1)

        self.frontend = self.context.socket(zmq.ROUTER)

        frontend_addr = 'tcp://*:{!s}'.format(frontend_port)

        self.frontend.bind(frontend_addr)

        logger.info('Binding frontend to %r', frontend_addr)

        self.backend = self.context.socket(zmq.ROUTER)

        backend_addr = 'tcp://*:{!s}'.format(backend_port)

        self.backend.bind(backend_addr)

        logger.info('Binding backend to %r', backend_addr)

        self.poll_both = zmq.Poller()

        self.poll_both.register(self.frontend, zmq.POLLIN)

        self.poll_both.register(self.backend, zmq.POLLIN)

        self.redis = redis.Redis(host=redis_host, db=1)

        logger.info('Connected to redis db %r', self.redis)

    def update_labels(self, yaml_data, labels):
        try:
            yaml_data['metadata']['labels'].update(labels)
        except AttributeError:
            # Labels not a dict, is this really a possible case?
            yaml_data['metadata']['labels'] = labels
        except KeyError:
            # Labels does not exists
            yaml_data['metadata'] = {
                'labels': labels,
            }

    def request_resources(self, request, resource_uuid, labels):
        requested = {}

        try:
            for item in request:
                yaml_data = yaml.safe_load(item)

                self.update_labels(yaml_data, labels)

                kind = yaml_data['kind']

                if kind == 'Pod':
                    self.core.create_namespaced_pod(body=yaml_data, namespace=self.namespace)

                    key = '{!s}:{!s}:Pod'.format(resource_uuid, yaml_data['metadata']['name'])
                elif kind == 'Deployment':
                    self.apps.create_namespaced_deployment(body=yaml_data, namespace=self.namespace)

                    key = '{!s}:{!s}:Deployment'.format(resource_uuid, yaml_data['metadata']['name'])
                elif kind == 'Service':
                    self.core.create_namespaced_service(body=yaml_data, namespace=self.namespace)

                    key = '{!s}:{!s}:Service'.format(resource_uuid, yaml_data['metadata']['name'])
                elif kind == 'Ingress':
                    self.extensions.create_namespaced_ingress(body=yaml_data, namespace=self.namespace)

                    key = '{!s}:{!s}:Ingress'.format(resource_uuid, yaml_data['metadata']['name'])
                else:
                    raise Exception('Requested an unsupported resource')

                requested[key] = yaml_data
        except client.rest.ApiException as e:
            if e.status == 409:
                logger.info('Resources already exists on the server')

                pass
            else:
                logger.error('Error allocating resources')

                raise ResourceAllocationError(yaml_data['metadata']['name'])

        return requested

    def wait_resources(self, resources, namespace, wait_timeout):
        for key, value in resources.items():
            logger.info('Waiting on resource %r', key)

            start = time.time()

            while True:
                kind = value['kind']

                if kind == 'Pod':
                    output = self.core.read_namespaced_pod_status(value['metadata']['name'], namespace)

                    if all([x.ready for x in output.status.container_statuses]):
                        logger.info('Pod %r containers are ready', value['metadata']['name'])

                        break
                elif kind == 'Deployment':
                    output = self.apps.read_namespaced_deployment_status(value['metadata']['name'], namespace)

                    if output.spec.replicas == output.status.ready_replicas:
                        logger.info('Deployment %r Pods are ready', value['metadata']['name'])

                        break
                else:
                    break

                if (time.time() - start) < wait_timeout:
                    raise ResourceTimeoutError(value['metadata']['name'])

                time.sleep(2)

    def allocate_resources(self, request):
        """ Allocate resources.

        Only support the following resources types: Pod, Deployment, Service, Ingress.

        Args:
            request: A list of YAML strs to create in the cluster.
        """
        logger.info('Allocating resources')

        resource_uuid = str(uuid4())[:8]

        # Create label with a uuid which can be queried for by the kube-monitor
        labels = {
            'app.kubernetes.io/resource-group-uuid': resource_uuid,
        }

        # Try to create the resources in the cluster, dont check to see if they
        # succeed, this might need to change and block until everything is up and
        # running
        try:
            requested = self.request_resources(request, resource_uuid, labels)
        except ResourceAllocationError as e:
            logger.error('Error allocation resource %s', e.name)

            raise e

        if self.wait_for_resources:
            try:
                self.wait_resources(requested, self.namespace, self.resource_timeout)
            except ResourceTimeoutError as e:
                logger.error('Timed out waiting for resource')

                raise e

        expired = (datetime.datetime.now() + datetime.timedelta(seconds=self.lifetime)).timestamp()

        for key in requested.keys():
            self.redis.hset('resource', key, expired)

    def handle_backend_frames(self, frames):
        address = frames[0]

        logger.debug('Handling frames from backend %r', frames[:3])

        if frames[1] == constants.RESOURCE:
            try:
                self.allocate_resources(json.loads(frames[2]))
            except (ResourceAllocationError, ResourceTimeoutError) as e:
                # Notify error in allocating resources
                new_frames = [address, constants.ERR, str(e)]
            else:
                # Allocate resources
                new_frames = [address, constants.ACK]

            self.backend.send_multipart(new_frames)
        elif frames[1] == constants.ACK:
            logger.info('Received ack from backend %r', address)

            self.waiting_ack.pop(address, None)
        else:
            logger.debug('Received message from backend %r', address)

            version = frames[2]

            self.workers.ready(Worker(address, version))

            if frames[1] == constants.READY:
                pass
            elif frames[1] == constants.HEARTBEAT:
                pass
            else:
                logger.error('Unknown message %r', frames)

    def handle_frontend_frames(self, frames):
        self.frontend.send_multipart(frames)

        version = frames[2]

        with self.redis.lock('job_queue'):
            self.redis.rpush(version, json_encoder(frames))

        logger.info('Added frames to %r queue %r', version, frames)

    def handle_unacknowledge_requests(self):
        for address, waiting in list(self.waiting_ack.items()):
            if time.time() >= waiting.expiry:
                logger.info('Backend %r never acknowledged request', address)

                self.waiting_ack.pop(address, None)

                # Reinsert the request frames infront of the queue
                with self.redis.lock('job_queue'):
                    self.redis.lpush(waiting.version, json_encoder(waiting.frames))

    def dispatch_workers(self):
        try:
            with self.redis.lock('job_queue', blocking_timeout=4):
                for address, worker in list(self.workers.queue.items()):
                    logger.debug('Processing waiting worker %r', address)

                    try:
                        frames = json_decoder(self.redis.lpop(worker.version))
                    except Exception:
                        logger.debug('No work found for version %r', worker.version)

                        continue

                    frames.insert(0, address)

                    frames.insert(1, constants.REQUEST)

                    logger.debug('Sending frames to worker %r', frames)

                    self.backend.send_multipart(frames)

                    self.waiting_ack[address] = WaitingAck(worker.version, frames[2:])

                    # Remove the worker
                    self.workers.queue.pop(address, None)

                    logger.info('Dispatch work to %r', address)
        except redis.lock.LockError:
            pass

    def monitor(self, frontend_port, backend_port, redis_host):
        """ Main loop for the load balancer.
        """
        self.initialize(frontend_port, backend_port, redis_host)

        while True:
            socks = dict(self.poll_both.poll(constants.HEARTBEAT_INTERVAL * 1000))

            if socks.get(self.backend) == zmq.POLLIN:
                frames = self.backend.recv_multipart()

                if not frames:
                    break

                self.handle_backend_frames(frames)

            if socks.get(self.frontend) == zmq.POLLIN:
                frames = self.frontend.recv_multipart()

                if not frames:
                    break

                self.handle_frontend_frames(frames)

            if self.workers.heartbeat_time:
                self.workers.heartbeat(self.backend)

            self.workers.purge()

            self.dispatch_workers()

            self.handle_unacknowledge_requests()


def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--log-level', help='Logging level', choices=logging._nameToLevel.keys(), default='INFO')

    parser.add_argument('--redis-host', help='Redis host', required=True)

    parser.add_argument('--frontend-port', help='Frontend port', type=int, default=7777)

    parser.add_argument('--backend-port', help='Backend port', type=int, default=7778)

    parser.add_argument('--namespace', help='Kubernetes namespace to monitor', default='default')

    parser.add_argument('--timeout', help='Resource monitor timeout', type=int, default=30)

    parser.add_argument('--resource-timeout', help='Time in seconds allowed for resources to be allocated', type=int,
                        default=240)

    parser.add_argument('--lifetime', help='Time in seconds to let resources live', type=int, default=3600)

    parser.add_argument('--wait-for-resources', help='The provisioner will wait until resources are allocated '
                        'before replying to the backend', type=bool, default=False)

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)

    provisioner = Provisioner(**vars(args))

    provisioner.start()

    monitor = kube_cluster.KubeCluster(args.redis_host, args.namespace, args.timeout, False)

    monitor.start()

    provisioner.join()

    monitor.join()


if __name__ == '__main__':
    main()