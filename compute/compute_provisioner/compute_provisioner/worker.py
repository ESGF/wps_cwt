import logging
import time

import zmq

from compute_provisioner import constants
from compute_provisioner import metrics

logger = logging.getLogger('compute_provisioner.worker')

# Handler types
REQUEST_TYPE = 'execute'
RESOURCE_TYPE = 'resource'
ERROR_TYPE = 'error'


class Worker(object):
    def __init__(self, version):
        self.version = version

        self.context = None

        self.worker = None

        self.poller = None

        self.heartbeat_at = None

        self.liveness = constants.HEARTBEAT_LIVENESS

        self.interval = constants.INTERVAL_INIT

        self.metrics = metrics.Metrics()

    def initialize(self, worker_addr):
        """ Initializes the worker.
        """
        self.worker = self.context.socket(zmq.DEALER)

        self.worker.connect('tcp://{!s}'.format(worker_addr))

        logger.info('Created dealer socket and connected to %r', worker_addr)

        self.poller.register(self.worker, zmq.POLLIN)

        self.worker.send_multipart([constants.READY, self.version])

        logger.info('Notifying queue, ready for work')

    def run(self, worker_addr, callback_handler):
        self.context = zmq.Context(1)

        self.poller = zmq.Poller()

        self.initialize(worker_addr)

        self.heartbeat_at = time.time() + constants.HEARTBEAT_INTERVAL

        pending_request = None

        while True:
            socks = dict(self.poller.poll(constants.HEARTBEAT_INTERVAL * 1000))

            if socks.get(self.worker) == zmq.POLLIN:
                frames = self.worker.recv_multipart()

                if not frames:
                    logger.error('Received empty message from queue')

                    self.metrics.inc('recv_queue_empty')

                    break

                logger.debug('Handling frames from provisioner %r', frames[:3])

                if len(frames) == 1 and frames[0] == constants.HEARTBEAT:
                    self.liveness = constants.HEARTBEAT_LIVENESS

                    logger.debug('Received heartbeat from queue, setting liveness to %r', self.liveness)

                    self.metrics.inc('recv_heartbeat')
                elif frames[0] == constants.REQUEST:
                    logger.info('Received request from queue %r', frames)

                    self.metrics.inc('recv_request')

                    pending_request = frames[3:]

                    resources = callback_handler(RESOURCE_TYPE, frames[3:])

                    logger.info('Resources %r', resources)

                    self.worker.send_multipart([constants.RESOURCE, resources.encode()])

                    self.metrics.inc('sent_resource')
                elif frames[0] == constants.ACK:
                    # Resource allocation ack, moving forward with execution
                    logger.info('Received data from queue %r', frames)

                    self.metrics.inc('recv_ack')

                    if pending_request is None:
                        logger.error('Recieved ack from resource allocation but have no pending request')
                    else:
                        callback_handler(REQUEST_TYPE, pending_request)

                        # Send provisioner ack
                        self.worker.send_multipart([constants.ACK])

                        self.metrics.inc('sent_ack')

                        # Clear out pending request
                        pending_request = None
                elif frames[0] == constants.ERR:
                    # Resource allocation error
                    logger.info('Received error from queue %r', frames[1])

                    callback_handler(ERROR_TYPE, pending_request + frames)

                    self.metrics.inc('recv_err')
                else:
                    logger.error('Received unknown message from queue %r', frames)

                    self.metrics.inc('recv_unknown')

                self.interval = constants.INTERVAL_INIT
            else:
                self.liveness -= 1

                if self.liveness == 0:
                    logger.error('Heartbeat failed cannot reach queue')
                    logger.error('Reconnect in %r', self.interval)

                    time.sleep(self.interval)

                    if self.interval < constants.INTERVAL_MAX:
                        self.interval *= 2

                    self.poller.unregister(self.worker)

                    self.worker.setsockopt(zmq.LINGER, 0)

                    self.worker.close()

                    self.initialize(worker_addr)

                    self.liveness = constants.HEARTBEAT_LIVENESS

                    self.metrics.inc('missed_heartbeat')

            if time.time() > self.heartbeat_at:
                self.heartbeat_at = time.time() + constants.HEARTBEAT_INTERVAL

                logger.debug('Sending heartbeat to queue')

                self.worker.send_multipart([constants.HEARTBEAT, self.version])

                self.metrics.inc('sent_heartbeat')


def demo():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--log-level', help='Logging level', choices=logging._nameToLevel.keys(), default='INFO')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)

    def handler(frames):
        logger.info('Received frames %r', frames)

    worker = Worker(b'devel')

    worker.run('0.0.0.0:7778', handler)
