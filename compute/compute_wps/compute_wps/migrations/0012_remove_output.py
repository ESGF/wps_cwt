# Generated by Django 3.1.2 on 2020-10-29 04:35

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('compute_wps', '0011_remove_process_backend'),
    ]

    operations = [
        migrations.DeleteModel(
            name='Output',
        ),
    ]