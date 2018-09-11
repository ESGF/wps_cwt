#! /usr/bin/env python

import json
import os
import shutil

import cwt
from celery.utils.log import get_task_logger
from django.conf import settings

from wps import metrics
from wps import models
from wps import WPSError
from wps.tasks import base

logger = get_task_logger('wps.tasks.job')

@base.cwt_shared_task()
def job_started(self, job_id):
    job = self.load_job(job_id)

    job.started()

@base.cwt_shared_task()
def job_succeeded(self, attrs, variables, output_path, move_path, var_name,
                  process_id, user_id, job_id):
    job = self.load_job(job_id)

    user = self.load_user(user_id)

    process = self.load_process(process_id)

    if move_path is not None:
        try:
            os.makedirs(os.path.dirname(move_path))
        except OSError:
            raise WPSError('Failed to create output directory')

        shutil.move(output_path, move_path)
        
        output_path = move_path

    relpath = os.path.relpath(output_path, settings.WPS_PUBLIC_PATH)

    url = settings.WPS_DAP_URL.format(filename=relpath)

    output = cwt.Variable(url, var_name)

    job.succeeded(json.dumps(output.parameterize()))

    process.track(user)

    for var in variables:
        models.File.track(user, var)

    return attrs
