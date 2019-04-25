#! /usr/bin/env python

import json

from celery.utils.log import get_task_logger
from django.conf import settings
from django.core.mail import EmailMessage

from wps import metrics
from wps import models
from wps.tasks import base

logger = get_task_logger('wps.tasks.job')

JOB_SUCCESS_MSG = """
Hello {name},

<pre>
Job {job.pk} has finished successfully.
</pre>

<pre>
{outputs}
</pre>

<pre>
<a href="{settings.WPS_JOB_URL}">View job history</a>
</pre>

<pre>
Thank you,
ESGF Compute Team
</pre>
"""

JOB_FAILED_MSG = """
Hello {name},

<pre>
Job {job.pk} has failed.
</pre>

<pre>
{error}
</pre>

Report errors on <a href="https://github.com/ESGF/esgf-compute-wps/issues/new?template=Bug_report.md">GitHub</a>.

<pre>
Thank you,
ESGF Compute Team
</pre>
"""


def send_failed_email(context, error):
    if context.user.first_name is None:
        name = context.user.username
    else:
        name = context.user.get_full_name()

    msg = JOB_FAILED_MSG.format(name=name, job=context.job, error=error)

    email = EmailMessage('Job Failed', msg, to=[context.user.email, ])

    email.content_subtype = 'html'

    email.send(fail_silently=True)


def send_success_email(context, variable):
    if context.user.first_name is None:
        name = context.user.username
    else:
        name = context.user.get_full_name()

    if not isinstance(variable, (list, tuple)):
        variable = [variable, ]

    outputs = '\n'.join('<a href="{!s}.html">{!s}|{!s}</a>'.format(
        x.uri,
        x.var_name,
        x.name)
        for x in variable)

    msg = JOB_SUCCESS_MSG.format(name=name, job=context.job, settings=settings, outputs=outputs)

    email = EmailMessage('Job Success', msg, to=[context.user.email, ])

    email.content_subtype = 'html'

    email.send(fail_silently=True)


def send_success_email_data(context, outputs):
    if context.user.first_name is None:
        name = context.user.username
    else:
        name = context.user.get_full_name()

    msg = JOB_SUCCESS_MSG.format(name=name, job=context.job, settings=settings, outputs=outputs)

    email = EmailMessage('Job Success', msg, to=[context.user.email, ])

    email.content_subtype = 'html'

    email.send(fail_silently=True)


@base.cwt_shared_task()
def job_started(self, context):
    context.job.started()

    return context


@base.cwt_shared_task()
def job_succeeded(self, context):
    if isinstance(context.output, str):
        context.job.succeeded(context.output)

        send_success_email_data(context, context.output)
    else:
        if len(context.output) == 1:
            context.job.succeeded(json.dumps(context.output[0]))
        else:
            context.job.succeeded(json.dumps(context.output))

        send_success_email(context, context.output)

    context.process.track(context.user)

    for input in context.inputs:
        models.File.track(context.user, input)

        metrics.track_file(input)

    return context
