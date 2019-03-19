from rest_framework import serializers
from rest_framework.reverse import reverse

from wps import models

class MessageSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Message
        fields = ('id', 'created_date', 'message', 'percent')

class StatusSerializer(serializers.ModelSerializer):
    messages = MessageSerializer(many=True)

    class Meta:
        model = models.Status
        fields = ('id', 'status', 'created_date', 'messages', 
                  'output', 'exception')

class StatusHyperlink(serializers.HyperlinkedRelatedField):
    view_name = 'status-detail'
    queryset = models.Status.objects.all()

    def get_url(self, obj, view_name, request, format):
        # Construct the url from the related job pk and the status pk
        url_kwargs = {
            'job_pk': obj.job.pk,
            'pk': obj.pk
        }

        return reverse(view_name, kwargs=url_kwargs, request=request, format=format)

class JobSerializer(serializers.ModelSerializer):
    server = serializers.SlugRelatedField(read_only=True,
                                          slug_field='host')

    process = serializers.SlugRelatedField(read_only=True,
                                           slug_field='identifier')

    elapsed = serializers.ReadOnlyField()

    latest_status = serializers.ReadOnlyField()

    accepted_on = serializers.ReadOnlyField()

    status = StatusHyperlink(many=True)

    class Meta:
        model = models.Job
        fields = ('id', 'server', 'process', 'extra', 'elapsed',
                  'latest_status', 'accepted_on', 'status')