# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-04-27 16:00
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('wps', '0011_remove_unused_oauth2_fields'),
    ]

    operations = [
        migrations.CreateModel(
            name='MPC',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('openid', models.TextField()),
                ('api_key', models.CharField(max_length=256)),
                ('cert', models.TextField()),
                ('user', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]