from __future__ import unicode_literals

import base64
import datetime
import hashlib
import json
import logging
import os
import random
import string
import time
from urlparse import urlparse

import cdms2
import cwt
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.db import transaction
from django.db.models import F
from django.db.models import signals
from django.db.models.query_utils import Q
from django.utils import timezone
from openid import association
from openid.store import interface
from openid.store import nonce

import wps
from wps import helpers
from wps import metrics

logger = logging.getLogger('wps.models')

ProcessAccepted = 'ProcessAccepted'
ProcessStarted = 'ProcessStarted'
ProcessPaused = 'ProcessPaused'
ProcessSucceeded = 'ProcessSucceeded'
ProcessFailed = 'ProcessFailed'

class DjangoOpenIDStore(interface.OpenIDStore):
    # Heavily borrowed from http://bazaar.launchpad.net/~ubuntuone-pqm-team/django-openid-auth/trunk/view/head:/django_openid_auth/store.py

    def storeAssociation(self, server_url, association):
        try:
            assoc = OpenIDAssociation.objects.get(server_url=server_url,
                                                 handle=association.handle)
        except OpenIDAssociation.DoesNotExist:
            assoc = OpenIDAssociation(server_url=server_url,
                                      handle=association.handle,
                                      secret=base64.encodestring(association.secret),
                                      issued=association.issued,
                                      lifetime=association.lifetime,
                                      assoc_type=association.assoc_type)
        else:
            assoc.secret = base64.encodestring(association.secret)
            assoc.issued = association.issued
            assoc.lifetime = association.lifetime
            assoc.assoc_type = association.assoc_type

        assoc.save()

    def getAssociation(self, server_url, handle=None):
        if handle is not None:
            assocs = OpenIDAssociation.objects.filter(server_url=server_url,
                                                     handle=handle)
        else:
            assocs = OpenIDAssociation.objects.filter(server_url=server_url)

        active = []
        expired = []

        for a in assocs:
            assoc = association.Association(a.handle,
                                            base64.decodestring(a.secret.encode('utf-8')),
                                            a.issued,
                                            a.lifetime,
                                            a.assoc_type)

            expires_in = assoc.getExpiresIn()

            if expires_in == 0:
                expired.append(a)
            else:
                active.append((assoc.issued, assoc))

        for e in expired:
            e.delete()

        if len(active) == 0:
            return None

        active.sort()

        return active[-1][1]

    def removeAssociation(self, server_url, handle):
        assocs = OpenIDAssociations.objects.filter(server_url=server_url, handle=handle)

        for a in assocs:
            a.delete()

        return len(assocs) > 0

    def useNonce(self, server_url, timestamp, salt):
        if abs(timestamp - time.time()) > nonce.SKEW:
            return False

        try:
            ononce = OpenIDNonce.objects.get(server_url__exact=server_url,
                                             timestamp__exact=timestamp,
                                             salt__exact=salt)
        except OpenIDNonce.DoesNotExist:
            ononce = OpenIDNonce.objects.create(server_url=server_url,
                                                timestamp=timestamp,
                                                salt=salt)
            return True

        return False

    def cleanupNonces(self):
        # Will implement later since it's not called
        pass

    def cleanupAssociations(self):
        # Will implement later since it's not called
        pass

class OpenIDNonce(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    server_url = models.CharField(max_length=2048)
    timestamp = models.IntegerField()
    salt = models.CharField(max_length=40)

    class Meta:
        unique_together = ('server_url', 'timestamp', 'salt')

class OpenIDAssociation(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    server_url = models.CharField(max_length=2048)
    handle = models.CharField(max_length=256)
    secret = models.TextField(max_length=256)
    issued = models.IntegerField()
    lifetime = models.IntegerField()
    assoc_type = models.TextField(max_length=64)

    class Meta:
        unique_together = ('server_url', 'handle')

class File(models.Model):
    name = models.CharField(max_length=256)
    host = models.CharField(max_length=256)
    variable = models.CharField(max_length=64)
    url = models.TextField()
    requested = models.PositiveIntegerField(default=0)

    class Meta:
        unique_together = ('name', 'host')

    @staticmethod
    def track(user, variable):
        url_parts = urlparse(variable.uri)

        parts = url_parts[2].split('/')

        file_obj, _ = File.objects.get_or_create(
            name=parts[-1],
            host=url_parts[1],
            variable=variable.var_name,
            url=variable.uri
        )

        file_obj.requested = F('requested') + 1

        file_obj.save()

        user_file_obj, _ = UserFile.objects.get_or_create(user=user, file=file_obj)

        user_file_obj.requested = F('requested') + 1

        user_file_obj.save()

    def to_json(self):
        return {
            'name': self.name,
            'host': self.host,
            'variable': self.variable,
            'url': self.url,
            'requested': self.requested
        }

class UserFile(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    file = models.ForeignKey(File, on_delete=models.CASCADE)
    requested_date = models.DateTimeField(auto_now=True)
    requested = models.PositiveIntegerField(default=0)

    def to_json(self):
        data = self.file.to_json()

        data['requested'] = self.requested

        data['requested_date'] = self.requested_date

        return data

class Cache(models.Model):
    uid = models.CharField(max_length=256)
    url = models.CharField(max_length=513)
    dimensions = models.TextField()
    added_date = models.DateTimeField(auto_now_add=True)
    accessed_date = models.DateTimeField(null=True)
    size = models.DecimalField(null=True, max_digits=16, decimal_places=8)

    @property
    def local_path(self):
        dimension_hash = hashlib.sha256(self.uid+self.dimensions).hexdigest()

        file_name = '{}.nc'.format(dimension_hash)

        return os.path.join(settings.WPS_CACHE_PATH, file_name)

    @property
    def valid(self):
        logger.info('Validating cached file %r', self.local_path)

        if not os.path.exists(self.local_path):
            logger.info('Invalid cached file %r is missing', self.local_path)

            return False

        try:
            data = helpers.decoder(self.dimensions)
        except ValueError:
            logger.info('Invalid failed to load dimensions of cached file %r', self.dimensions)

            return False

        try:
            var_name = data['var_name']
        except KeyError:
            logger.info('Invalid missing var_name from cached metadata')

            return False
       
        del data['var_name']

        infile = None

        try:
            infile = cdms2.open(self.local_path)
        except cdms2.CDMSError:
            logger.info('Invalid failed to open cached file %r', self.local_path)

            return False
        else:
            if var_name not in infile.variables:
                logger.info('Invalid variable %r not in cached file', var_name)

                return False

            variable = infile[var_name]

            for name, value in data.iteritems():
                axis_index = variable.getAxisIndex(name)

                if axis_index == -1:
                    logger.info('Invalid axis %r is missing from cached file', name)

                    return False

                axis = variable.getAxis(axis_index)

                expected = value.stop - value.start

                if expected != axis.shape[0]:
                    logger.info('Invalid axis %r expected %r found %r', axis.id, expected, axis.shape[0])

                    return False
        finally:
            if infile:
                infile.close()

        return True

    def is_superset(self, domain):
        try:
            cached = json.loads(self.dimensions, object_hook=helpers.json_loads_object_hook)
        except ValueError:
            return False

        del cached['var_name']

        for name, value in domain.iteritems():
            try:
                cached_value = cached[name] 
            except KeyError:
                return False

            if (value.start < cached_value.start or
                    value.stop > cached_value.stop):
                return False

        return True

    def accessed(self, when=None):
        if when is None:
            when = timezone.now()

        self.accessed = when

        self.save()

    def set_size(self):
        stat = os.stat(self.local_path)

        self.size = stat.st_size

        self.save()

    def estimate_size(self):
        with cdms2.open(self.url) as infile:
            var_name, dimensions = self.dimensions.split('!')

            size = reduce(lambda x, y: x * y, infile[var_name].shape)

            size = size * infile[var_name].dtype.itemsize

        return size

class Auth(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)

    openid = models.TextField()
    openid_url = models.CharField(max_length=256)
    type = models.CharField(max_length=64)
    cert = models.TextField()
    api_key = models.CharField(max_length=128)
    extra = models.TextField()

    def update(self, auth_type, certs, **kwargs):
        if self.api_key == '':
            self.api_key = ''.join(random.choice(string.ascii_letters+string.digits) for _ in xrange(64))

        self.type = auth_type

        self.cert = ''.join(certs)

        extra = json.loads(self.extra or '{}')

        extra.update(kwargs)

        self.extra = json.dumps(extra)

        self.save()

class Process(models.Model):
    identifier = models.CharField(max_length=128, unique=True)
    backend = models.CharField(max_length=128)
    abstract = models.TextField()
    description = models.TextField()
    enabled = models.BooleanField(default=True)

    def track(self, user):
        user_process_obj, _ = UserProcess.objects.get_or_create(user=user, process=self)

        user_process_obj.requested = F('requested') + 1

        user_process_obj.save()

    def to_json(self, usage=False):
        return {
            'identifier': self.identifier,
            'backend': self.backend,
            'enabled': self.enabled
        }

class UserProcess(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    process = models.ForeignKey(Process, on_delete=models.CASCADE)
    requested_date = models.DateTimeField(auto_now=True)
    requested = models.PositiveIntegerField(default=0)

    def to_json(self):
        data = self.process.to_json()

        data['requested'] = self.requested

        data['requested_date'] = self.requested_date

        return data

class Server(models.Model):
    host = models.CharField(max_length=128, unique=True)
    added_date = models.DateTimeField(auto_now_add=True)
    status = models.IntegerField(default=1)
    capabilities = models.TextField()

    processes = models.ManyToManyField(Process)

class Job(models.Model):
    server = models.ForeignKey(Server, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    process = models.ForeignKey(Process, on_delete=models.CASCADE, null=True)
    steps_completed = models.IntegerField(default=0)
    steps_total = models.IntegerField(default=0)
    extra = models.TextField(null=True)

    @property
    def details(self):
        return {
            'id': self.id,
            'server': self.server.host,
            'process': self.process.identifier,
            'elapsed': self.elapsed,
            #'status': self.status
        }

    @property
    def status(self):
        return [
            {
                'created_date': x.created_date,
                'updated_date': x.updated_date,
                'status': x.status,
                'exception': x.exception,
                'output': x.output,
                'messages': [y.details for y in x.message_set.all().order_by('created_date')]
            } for x in self.status_set.all().order_by('created_date')
        ]

    @property
    def elapsed(self):
        started = self.status_set.filter(status='ProcessStarted')

        ended = self.status_set.filter(Q(status='ProcessSucceeded') | Q(status='ProcessFailed'))

        if len(started) == 0 or len(ended) == 0:
            return 'No elapsed'

        elapsed = ended[0].created_date - started[0].created_date

        return '{}.{}'.format(elapsed.seconds, elapsed.microseconds)

    @property
    def report(self):
        location = settings.WPS_STATUS_LOCATION.format(job_id=self.id)

        latest = self.status_set.latest('created_date')

        if latest.exception is not None:
            ex_report = cwt.wps.CreateFromDocument(latest.exception)

            status = cwt.wps.status_failed_from_report(ex_report)
        elif latest.status == ProcessAccepted:
            status = cwt.wps.status_accepted('Job accepted')
        elif latest.status == ProcessStarted:
            message = latest.message_set.latest('created_date')

            status = cwt.wps.status_started(message.message, message.percent)
        elif latest.status == ProcessPaused:
            message = latest.message_set.latest('created_date')

            status = cwt.wps.status_paused(message.message, message.percent)
        elif latest.status == ProcessSucceeded:
            status = cwt.wps.status_succeeded('Job succeeded')

        outputs = []

        outputs.append(cwt.wps.output_data('output', 'output', latest.output))

        process = cwt.wps.process(self.process.identifier, self.process.identifier, '1.0.0')

        args = [
            process,
            status,
            '1.0.0',
            'en-US',
            settings.WPS_ENDPOINT,
            location,
            outputs
        ]

        response = cwt.wps.execute_response(*args)

        cwt.bds.reset()

        return response.toxml(bds=cwt.bds)

    @property
    def is_started(self):
        latest = self.status_set.latest('created_date')

        return latest is not None and latest.status == ProcessStarted

    @property
    def steps_progress(self):
        try:
            return (self.steps_completed + 1) * 100.0 / self.steps_total
        except ZeroDivisionError:
            return 0.0

    @property
    def steps(self):
        return self.steps_total

    @steps.setter
    def steps(self, x):
        self.steps_total = x

        self.save()

    def steps_inc_total(self, steps=None):
        if steps is None:
            steps = 1

        self.steps_total = F('steps_total') + steps

        self.save()

    def steps_reset(self):
        self.steps_total = 0

        self.steps_completed = 0

        self.save()

    def steps_inc(self, steps=None):
        if steps is None:
            steps = 1

        self.steps_completed = F('steps_completed') + steps

        self.save()

    def accepted(self):
        self.status_set.create(status=ProcessAccepted)

        metrics.JOBS_IN_QUEUE.inc()

        metrics.JOBS_QUEUED.labels(self.process.identifier).inc()

    def started(self):
        if not self.is_started:
            status = self.status_set.create(status=ProcessStarted)

            status.set_message('Job Started', self.steps_progress)

            metrics.JOBS_IN_QUEUE.dec()

            metrics.JOBS_RUNNING.inc()

    def succeeded(self, output=None):
        if output is None:
            output = ''

        status = self.status_set.create(status=ProcessSucceeded)

        status.output = output

        status.save()

        metrics.JOBS_RUNNING.dec()

        metrics.JOBS_COMPLETED.labels(self.process.identifier).inc()

    def failed(self, exception=None):
        if self.is_started:
            self.update('Job Failed')

        if exception is None:
            exception = 'Unknown reason'

        with transaction.atomic():
            status = self.status_set.create(status=ProcessFailed)

            status.exception = wps.exception_report(exception, cwt.ows.NoApplicableCode)

            status.save()

            metrics.JOBS_RUNNING.dec()

            metrics.JOBS_FAILED.labels(self.process.identifier).inc()

    def retry(self, exception):
        self.update('Retrying... {}', exception)

    def update(self, message, *args, **kwargs):
        message = message.format(*args, **kwargs)

        started = self.status_set.filter(status=ProcessStarted).latest('created_date')

        started.set_message(message, self.steps_progress)

    def statusSince(self, date):
        return [
            {
                'created_date': x.created_date,
                'updated_date': x.updated_date,
                'status': x.status,
                'exception': x.exception,
                'output': x.output,
                'messages': [y.details for y in x.message_set.all().order_by('created_date')]
            } for x in self.status_set.filter(updated_date__gt=date)
        ]

class Status(models.Model):
    job = models.ForeignKey(Job, on_delete=models.CASCADE)

    created_date = models.DateTimeField(auto_now_add=True)
    updated_date = models.DateTimeField(auto_now_add=True)
    status = models.CharField(max_length=128)
    exception = models.TextField(null=True)
    output = models.TextField(null=True)

    def set_message(self, message, percent=None):
        self.message_set.create(message=message, percent=percent)

        self.updated_date = timezone.now()

        self.save()

class Message(models.Model):
    status = models.ForeignKey(Status, on_delete=models.CASCADE)

    created_date = models.DateTimeField(auto_now_add=True)
    percent = models.PositiveIntegerField(null=True)
    message = models.TextField(null=True)

    @property
    def details(self):
        return {
            'created_date': self.created_date,
            'percent': self.percent or 0,
            'message': self.message
        }
