import argparse
import json
import logging
import os
from functools import partial

import cwt
import jinja2

from compute_tasks import base
from compute_tasks import cdat
from compute_tasks import celery_ as celery
from compute_tasks import WPSError
from compute_tasks.job import job_started
from compute_tasks.job import job_succeeded
from compute_tasks.context import state_mixin
from compute_provisioner.worker import Worker
from compute_provisioner.worker import REQUEST_TYPE
from compute_provisioner.worker import RESOURCE_TYPE
from compute_provisioner.worker import ERROR_TYPE

logger = logging.getLogger('compute_tasks.backend')

# Set INGRESS_QUEUE to prevent breaking old code
DEFAULT_QUEUE = {
    'queue': 'ingress',
    'exchange': 'ingress',
    'routing_key': 'ingress',
}

QUEUE = {
    'cdat': {
        'queue': 'ingress',
        'exchange': 'ingress',
        'routing_key': 'ingress',
    },
    'default': {
        'queue': 'default',
        'exchange': 'default',
        'routing_key': 'default',
    },
}

QUEUE.update({'ingress': DEFAULT_QUEUE})


def queue_from_identifier(identifier):
    module, name = identifier.split('.')

    return QUEUE.get(module.lower(), DEFAULT_QUEUE)


def fail_job(state, job, e):
    try:
        state.job = job

        state.failed(str(e))
    except Exception:
        pass
    finally:
        state.job = None


def format_frames(frames):
    logger.debug('Handling frames %r', frames)

    version, identifier, data_inputs, job, user, process = [x.decode() for x in frames]

    data_inputs = celery.decoder(data_inputs)

    logger.info('Building celery workflow')

    extra = {
        'DASK_SCHEDULER': 'dask-scheduler-{!s}.{!s}.svc:8786'.format(user, os.environ['NAMESPACE']),
    }

    return identifier, data_inputs, job, user, process, extra


def validate_process(process):
    v = cdat.VALIDATION[process.identifier]

    if len(process.inputs) < v['min'] or len(process.inputs) > v['max']:
        raise WPSError('Validation failed, expected the number of inputs to be between {!s} and {!s}', v['min'], v['max'])

    for key, param_v in v['params'].items():
        required = param_v['required']

        param = process.parameters.get(key, None)

        if param is None and required:
            raise WPSError('Validation failed, expected parameter {!r}', key)

        if param is not None:
            type = param_v['type']

            for value in param.values:
                try:
                    type(value)
                except ValueError:
                    raise WPSError('Validation failed, could not convert parameter {!r} value {!r} to type {!r}', key, value, type.__name__)



def build_workflow(frames):
    data_inputs = frames[1]

    for process in data_inputs['operation']:
        validate_process(cwt.Process.from_dict(json.loads(process)))

    started = job_started.s(*frames).set(**DEFAULT_QUEUE)

    logger.info('Created job started task %r', started)

    queue = queue_from_identifier(frames[0])

    logger.info('Using queue %r', queue)

    process = base.get_process(frames[0]).s().set(**queue)

    logger.info('Found process %r for %r', frames[4], frames[0])

    succeeded = job_succeeded.s().set(**DEFAULT_QUEUE)

    logger.info('Created job stopped task %r', succeeded)

    return started | process | succeeded


def request_handler(frames, state):
    try:
        frames = format_frames(frames)

        workflow = build_workflow(frames)

        logger.info('Built workflow %r', workflow)

        workflow.delay()

        logger.info('Executed workflow')
    except Exception as e:
        logger.exception('Error executing celery workflow %r', e)

        fail_job(state, frames[2], e)


def resource_request(frames, env):
    version, identifier, data_inputs, job, user, process = [x.decode() for x in frames]

    resources = []

    data = {
        'image': os.environ['IMAGE'],
        'image_pull_secret': os.environ.get('IMAGE_PULL_SECRET', None),
        'image_pull_policy': os.environ.get('IMAGE_PULL_POLICY', 'Always'),
        'scheduler_cpu': os.environ.get('SCHEDULER_CPU', 1),
        'scheduler_memory': os.environ.get('SCHEDULER_MEMORY', '1Gi'),
        'worker_cpu': os.environ.get('WORKER_CPU', 1),
        'worker_memory': os.environ.get('WORKER_MEMORY', '1Gi'),
        'worker_nthreads': os.environ.get('WORKER_NTHREADS', 4),
        'traffic_type': os.environ.get('TRAFFIC_TYPE', 'development'),
        'dev': os.environ.get('DEV', False),
        'user': user,
        'workers': os.environ['WORKERS'],
        'data_claim_name': os.environ.get('DATA_CLAIM_NAME', 'data-pvc'),
    }

    dask_scheduler_pod = env.get_template('dask-scheduler-pod.yaml')

    dask_scheduler_service = env.get_template('dask-scheduler-service.yaml')

    dask_scheduler_ingress = env.get_template('dask-scheduler-ingress.yaml')

    dask_worker_deployment = env.get_template('dask-worker-deployment.yaml')

    resources.append(dask_scheduler_pod.render(**data))

    resources.append(dask_scheduler_service.render(**data))

    resources.append(dask_scheduler_ingress.render(**data))

    resources.append(dask_worker_deployment.render(**data))

    return json.dumps(resources)


def error_handler(frames, state):
    version, identifier, data_inputs, job, user, process = [x.decode() for x in frames[:-2]]

    logger.error('Resource allocation failed %r', frames[-1])

    fail_job(state, job, frames[-1].decode())


def load_processes(state, register_tasks=True):
    if register_tasks:
        for item in base.discover_processes():
            try:
                state.register_process(**item)
            except state_mixin.ProcessExistsError:  # pragma: no cover
                logger.info('Process %r already exists', item['identifier'])

                pass

    base.build_process_bindings()


def register_processes():
    parser = argparse.ArgumentParser()

    parser.add_argument('--log-level', help='Logging level', choices=logging._nameToLevel.keys(), default='INFO')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)

    state = state_mixin.StateMixin()

    state.init_api()

    load_processes(state)


def backend_argparse():  # pragma: no cover
    parser = argparse.ArgumentParser()

    parser.add_argument('--log-level', help='Logging level', choices=logging._nameToLevel.keys(), default='INFO')

    parser.add_argument('--queue-host', help='Queue to communicate with')

    parser.add_argument('-d', help='Development mode', action='store_true')

    return parser


def main():
    parser = backend_argparse()

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)

    logger.debug('CLI Arguments %r', args)

    logger.info('Loading templates')

    env = jinja2.Environment(loader=jinja2.PackageLoader('compute_tasks', 'templates'))

    state = state_mixin.StateMixin()

    state.init_api()

    base.discover_processes()

    base.build_process_bindings()

    # Need to get the supported version from somewhere
    # environment variable or hard code?
    worker = Worker(b'devel')

    queue_host = args.queue_host or os.environ['PROVISIONER_BACKEND']

    request_handler_partial = partial(request_handler, state=state)

    request_resource_partial = partial(resource_request, env=env)

    error_handler_partial = partial(error_handler, state=state)

    def callback_handler(type, frames):  # pragma: no cover
        if type == REQUEST_TYPE:
            value = request_handler_partial(frames)
        elif type == RESOURCE_TYPE:
            value = request_resource_partial(frames)
        elif type == ERROR_TYPE:
            value = error_handler_partial(frames)
        else:
            logger.error('Could not handle unknown type %r', type)

            raise Exception('Could not handle unknown type {!r}'.format(type))

        return value

    worker.run(queue_host, callback_handler)
