# -*- coding: utf-8 -*-
# Generated by Django 1.11.7 on 2018-01-04 01:22
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('wps', '0034_process_abstract'),
    ]

    operations = [
        migrations.AlterField(
            model_name='cache',
            name='accessed_date',
            field=models.DateTimeField(null=True),
        ),
    ]
