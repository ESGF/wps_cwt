# -*- coding: utf-8 -*-
# Generated by Django 1.11.25 on 2020-10-19 22:18
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('compute_wps', '0002_job_expiration'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='openidassociation',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='openidassociation',
            name='user',
        ),
        migrations.AlterUniqueTogether(
            name='openidnonce',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='openidnonce',
            name='user',
        ),
        migrations.DeleteModel(
            name='OpenIDAssociation',
        ),
        migrations.DeleteModel(
            name='OpenIDNonce',
        ),
    ]