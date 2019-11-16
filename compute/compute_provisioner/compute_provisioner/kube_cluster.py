import logging
import threading
import time

import redis
from kubernetes import client
from kubernetes import config

logger = logging.getLogger('provisioner.kube_cluster')


class KubeCluster(threading.Thread):
    def __init__(self, redis_host, namespace, timeout, dry_run):
        super(KubeCluster, self).__init__(target=self.monitor)

        self.redis_host = redis_host

        self.namespace = namespace

        self.timeout = timeout

        self.dry_run = dry_run

        config.load_incluster_config()

        self.core = client.CoreV1Api()

        self.apps = client.AppsV1Api()

        self.exts = client.ExtensionsV1beta1Api()

        self.redis = redis.Redis(self.redis_host, db=1)

    def check_rogue_expired_resources(self, kind, selector, resource, delete_func):
        logger.info('Checking %r %r resources', len(resource.items), kind)

        for x in resource.items:
            key = '{!s}:{!s}:{!s}'.format(x.metadata.labels[selector], x.metadata.name, kind)

            expire = self.redis.hget('resource', key)

            if expire is None:
                logger.info('Removing rogue resource %r %r', kind, x.metadata.name)

                if not self.dry_run:
                    delete_func(x.metadata.name, self.namespace)

                    logger.info('Removed %r %r', kind, x.metadata.name)
            else:
                expire = float(expire.decode())

                if expire < time.time():
                    logger.info('Removing expired resource %r %r', kind, x.metadata.name)

                    if not self.dry_run:
                        delete_func(x.metadata.name, self.namespace)

                        logger.info('Removed %r %r', kind, x.metadata.name)

                        self.redis.hdel('resource', key)
                else:
                    logger.debug('Found validate resource %r %r', kind, x.metadata.name)

    def check_resources(self):
        for key, expire in self.redis.hscan_iter('resource'):
            logger.debug('Checking key %r expire %r', key, expire)

            uuid, name, kind = key.decode().split(':')

            selector = 'app.kubernetes.io/resource-group-uuid={!s}'.format(uuid)

            if kind == 'Pod':
                resource = self.core.list_namespaced_pod(self.namespace, label_selector=selector)

                delete_func = self.core.delete_namespaced_pod
            elif kind == 'Deployment':
                resource = self.apps.list_namespaced_deployment(self.namespace, label_selector=selector)

                delete_func = self.apps.delete_namespaced_deployment
            elif kind == 'Service':
                resource = self.core.list_namespaced_service(self.namespace, label_selector=selector)

                delete_func = self.core.delete_namespaced_service
            elif kind == 'Ingress':
                resource = self.exts.list_namespaced_ingress(self.namespace, label_selector=selector)

                delete_func = self.exts.delete_namespaced_ingress
            else:
                logger.error('Cannot handle resource %r kind', kind)

                continue

            expire = float(expire)

            if len(resource.items) == 0:
                logger.info('Resource not found, removing key %r', name)

                # Resource not found
                if not self.dry_run:
                    self.redis.hdel('resource', key)
            elif time.time() > expire:
                logger.info('Resource expired, removing %r %r', kind, name)

                # Resource expired
                if not self.dry_run:
                    delete_func(name, self.namespace)

                    self.redis.hdel('resource', key)

        logger.info('Done checking for resources')

    def validate_existing_resources(self):
        selector = 'app.kubernetes.io/resource-group-uuid'

        pods = self.core.list_namespaced_pod(self.namespace, label_selector=selector)

        self.check_rogue_expired_resources('Pod', selector, pods, self.core.delete_namespaced_pod)

        deployments = self.apps.list_namespaced_deployment(self.namespace, label_selector=selector)

        self.check_rogue_expired_resources('Deployment', selector, deployments, self.apps.delete_namespaced_deployment)

        services = self.core.list_namespaced_service(self.namespace, label_selector=selector)

        self.check_rogue_expired_resources('Service', selector, services, self.core.delete_namespaced_service)

        ingresses = self.exts.list_namespaced_ingress(self.namespace, label_selector=selector)

        self.check_rogue_expired_resources('Ingress', selector, ingresses, self.exts.delete_namespaced_ingress)

    def monitor(self):
        self.validate_existing_resources()

        while True:
            self.check_resources()

            time.sleep(self.timeout)