import json
import os

import cdms2
import cwt
import requests
from cdms2 import MV2 as MV
from celery.task.control import inspect
from celery.utils.log import get_task_logger
from django.conf import settings
from django.db.models import Q
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


def query_single_value(type=int, **kwargs):
    try:
        data = query_prometheus(**kwargs)[0]
    except IndexError:
        return type()

    try:
        return type(data['value'][1])
    except (KeyError, IndexError):
        return type()

def query_multiple_value(key, type=int, **kwargs):
    results = {}

    data = query_prometheus(**kwargs)

    for item in data:
        try:
            name = item['metric'][key]
        except (KeyError, TypeError):
            continue

        try:
            value = item['value'][1]
        except (KeyError, IndexError):
            result[name] = type()
        else:
            results[name] = type(value)

    return results


@base.register_process('CDAT.metrics', abstract="""
                       Returns the current metrics of the server.
                       """, data_inputs=[], metadata={'inputs': 0})
@base.cwt_shared_task()
def metrics_task(self, user_id, job_id, **kwargs):
    job = self.load_job(job_id)

    user_jobs_queued = models.Job.objects.filter(status__status=models.ProcessAccepted).exclude(status__status=models.ProcessStarted).exclude(
        status__status=models.ProcessFailed).exclude(status__status=models.ProcessSucceeded).count()

    user_jobs_running = models.Job.objects.filter(status__status=models.ProcessStarted).exclude(
        status__status=models.ProcessFailed).exclude(status__status=models.ProcessSucceeded).count()

    operator_count = query_multiple_value('request', type=float, query='sum(wps_request_seconds_count) by (request)')

    operator_avg_time = query_multiple_value('request', type=float,
                                             query='avg(wps_request_seconds_sum) by (request)')

    operator = {}

    for item in set(operator_count.keys()+operator_avg_time.keys()):
        operator[item] = {}

        if item in operator_count:
            operator[item]['count'] = operator_count[item]

        if item in operator_avg_time:
            operator[item]['avg_time'] = operator_avg_time[item]

    file_count = query_multiple_value('url', query='sum(wps_file_accessed{url!=""}) by (url)')

    file = {}

    for item in file_count.keys():
        logger.info('%r', item)

        try:
            url_obj = models.File.objects.filter(url=item)[0]
        except IndexError:
            count = 0
        else:
            count = url_obj.userfile_set.all().distinct('user').count()

        file[item] = {'count': file_count[item], 'unique_users': count }

    data = {
        'health': {
            'jobs_running': query_single_value(type=int,
                                               query='sum(wps_jobs_running)'),
            'jobs_queued': query_single_value(type=int,
                                              query='sum(wps_jobs_in_queue)'),
            'user_jobs_running': user_jobs_running,
            'user_jobs_queued': user_jobs_queued,
            'cpu_avg': query_single_value(type=float,
                                          query='sum(rate(container_cpu_usage_seconds_total{namespace="default",container_name=~".*(ingress|wps).*"}[5m]))'),
            'cpu_count': query_single_value(type=int, query='sum(machine_cpu_cores)'),
            'memory_usage_avg_5m': query_single_value(type=float,
                                                      query='sum(avg_over_time(container_memory_usage_bytes{container_name=~".*(celery|wps).*"}[5m]))'),
            'memory_usage': query_single_value(type=float,
                                             query='sum(container_memory_usage_bytes{container_name=~".*(celery|wps).*"})'),
            'memory_available': query_single_value(type=int,
                                                   query='sum(container_memory_max_usage_bytes{container_name=~".*(celery|wps).*"})'),
            'wps_requests': query_single_value(type=int,
                                               query='sum(wps_request_seconds_count)'),
            'wps_requests_avg_5m': query_single_value(type=float,
                                                   query='sum(avg_over_time(wps_request_seconds_count[5m]))'),
        },
        'usage': {
            'files': file,
            'operators': operator,
            'output': query_single_value(type=float,
                                           query='sum(wps_process_bytes/wps_process_seconds >= 0)'),
            'local': query_single_value(type=float,
                                        query='sum(wps_cache_bytes/wps_cache_seconds >= 0)'),
            'download': query_single_value(type=float,
                                         query='sum(wps_ingress_bytes/wps_ingress_seconds >= 0)'),
        },
        'time': timezone.now().ctime(),
    }

    job.succeeded(json.dumps(data))

    return data
