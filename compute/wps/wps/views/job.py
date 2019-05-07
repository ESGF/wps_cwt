#! /usr/bin/env python

from urllib.parse import urlparse

from django import db
from django.db.models import F
from rest_framework import mixins
from rest_framework import viewsets
from rest_framework.authentication import BasicAuthentication
from rest_framework.response import Response
from rest_framework.permissions import DjangoModelPermissions

from wps import models
from wps import serializers


class InternalUserFileViewSet(mixins.CreateModelMixin,
                              viewsets.GenericViewSet):
    queryset = models.UserFile.objects.all()
    serializer_class = serializers.UserFileSerializer
    authentication_classes = (BasicAuthentication, )
    permission_classes = (DjangoModelPermissions, )

    def create(self, request, *args, **kwargs):
        try:
            user = models.User.objects.get(pk=kwargs['user_pk'])
        except models.User.DoesNotExist:
            return Response('User does not exist', status=400)

        parts = urlparse(request.data['url'])

        path_parts = parts.path.split('/')

        file, created = models.File.objects.get_or_create(name=path_parts[-1], host=parts.netloc)

        fields = ['requested', ]

        if created:
            file.url = request.data['url']

            file.variable = request.data['var_name']

            fields.extend(['url', 'variable'])

        file.requested = F('requested') + 1

        file.save(update_fields=fields)

        user_file, _ = models.UserFile.objects.get_or_create(user=user, file=file)

        user_file.requested = F('requested') + 1

        user_file.save(update_fields=['requested'])

        user_file.refresh_from_db()

        user_file_serializer = serializers.UserFileSerializer(instance=user_file)

        return Response(user_file_serializer.data, status=201)


class InternalUserProcessViewSet(mixins.CreateModelMixin,
                                 viewsets.GenericViewSet):
    queryset = models.UserProcess.objects.all()
    serializer_class = serializers.UserProcessSerializer
    authentication_classes = (BasicAuthentication, )
    permission_classes = (DjangoModelPermissions, )

    def create(self, request, *args, **kwargs):
        try:
            user = models.User.objects.get(pk=kwargs['user_pk'])
        except models.User.DoesNotExist:
            return Response('User does not exist', status=400)

        user_process, _ = models.UserProcess.objects.get_or_create(user=user, process=kwargs['process_pk'])

        user_process.requested = F('requested') + 1

        user_process.save(update_fields=['requested'])

        user_process.refresh_from_db()

        user_process_serializer = serializers.UserProcessSerializer(instance=user_process)

        return Response(user_process_serializer.data, status=201)


class InternalProcessViewSet(mixins.CreateModelMixin,
                             mixins.ListModelMixin,
                             viewsets.GenericViewSet):
    queryset = models.Process.objects.all()
    serializer_class = serializers.ProcessSerializer
    authentication_classes = (BasicAuthentication, )
    permission_classes = (DjangoModelPermissions, )


class InternalMessageViewSet(mixins.CreateModelMixin,
                             viewsets.GenericViewSet):
    queryset = models.Message.objects.all()
    serializer_class = serializers.MessageSerializer
    authentication_classes = (BasicAuthentication, )
    permission_classes = (DjangoModelPermissions, )

    def create(self, request, *args, **kwargs):
        try:
            status = models.Status.objects.get(job__pk=kwargs['job_pk'], pk=kwargs['status_pk'])
        except models.Status.DoesNotExist:
            return Response('Status does not exist', status=400)

        message_serializer = serializers.MessageSerializer(data=request.data)

        message_serializer.is_valid(raise_exception=True)

        message_serializer.save(status=status)

        return Response(message_serializer.data, status=201)


class InternalStatusViewSet(mixins.CreateModelMixin,
                            viewsets.GenericViewSet):
    queryset = models.Status.objects.all()
    serializer_class = serializers.StatusSerializer
    authentication_classes = (BasicAuthentication, )
    permission_classes = (DjangoModelPermissions, )

    def create(self, request, *args, **kwargs):
        try:
            job = models.Job.objects.get(pk=kwargs['job_pk'])
        except models.Job.DoesNotExist:
            return Response('Job does not exist', status=400)

        status_serializer = serializers.StatusSerializer(data=request.data)

        status_serializer.is_valid(raise_exception=True)

        try:
            status_serializer.save(job=job)
        except db.IntegrityError:
            return Response('Status already exists', status=400)

        return Response(status_serializer.data, status=201)


class StatusViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = models.Status.objects.all()
    serializer_class = serializers.StatusSerializer

    def get_queryset(self):
        user = self.request.user

        job_pk = self.kwargs['job_pk']

        return models.Status.objects.filter(job__pk=job_pk,
                                            job__user=user.pk)


class JobViewSet(mixins.ListModelMixin,
                 mixins.RetrieveModelMixin,
                 mixins.DestroyModelMixin,
                 viewsets.GenericViewSet):
    queryset = models.Job.objects.all()
    serializer_class = serializers.JobSerializer

    def get_queryset(self):
        user = self.request.user

        return models.Job.objects.filter(user=user.pk)
