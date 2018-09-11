import json
import os

import cdms2
import cwt
import requests
from cdms2 import MV2 as MV
from celery.task.control import inspect
from celery.utils.log import get_task_logger
from django.conf import settings
from django.utils import timezone

from wps import helpers
from wps import metrics
from wps import models
from wps import WPSError
from wps.tasks import base

logger = get_task_logger('wps.tasks.metrics')

def query_prometheus(**kwargs):
    response = requests.get(settings.METRICS_HOST, params=kwargs)

    if not response.ok:
        raise WPSError('Failed querying "{}" {}: {}', settings.METRICS_HOST,
                       response.reason, response.status_code)

    data = response.json()

    try:
        status = data['status']
    except KeyError:
        raise WPSError('Excepted JSON from prometheus request')

    logger.info('%r', data)

    return data['data']['result']

def query_single_value(type=str, **kwargs):
    try:
        data = query_prometheus(**kwargs)[0]
    except IndexError:
        return None

    return type(data['value'][1])

@base.register_process('CDAT.metrics', abstract="""
                       Returns the current metrics of the server.
                       """, data_inputs=[], metadata={'inputs': 0})
@base.cwt_shared_task()
def health(self, user_id, job_id, **kwargs):
    job = self.load_job(job_id)

    data = {
        'health': {
            'jobs_running': 0,
            'jobs_queued': 0,
            'users_running': 0,
            'users_queued': 0,
            'cpu_avg': 0,
            'cpu_min': 0,
            'cpu_max': 0,
            'cpu_count': query_single_value(type=int, query='sum(machine_cpu_cores)'),
            'memory_avg': 0,
            'memory_min': 0,
            'memory_max': 0,
            'memory_available': query_single_value(type=int,
                                                   query='sum(machine_memory_bytes)'),
            'disk_free_space': 0,
            'wps_requests_avg': 0,
        },
        'usage': {
            'files': {
            },
            'operators': {
            },
            'download': 0,
            'local': 0,
            'output': 0,
        },
        'time': timezone.now().ctime(),
    }

    job.succeeded(json.dumps(data))

    return data
