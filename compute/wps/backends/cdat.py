import json
import logging
import os
import re
import uuid

import celery
import cwt
import numpy as np
from django.conf import settings

from wps import helpers
from wps import models
from wps import tasks
from wps import WPSError
from wps.backends import backend
from wps.context import OperationContext
from wps.context import WorkflowOperationContext
from wps.tasks import base

logger = logging.getLogger('wps.backends')

class CDAT(backend.Backend):
    def populate_processes(self):
        logger.info('Registering processes for backend "local"')

        for name, proc in base.REGISTRY.iteritems():
            self.add_process(name, name.title(), metadata=proc.METADATA,
                             data_inputs=proc.DATA_INPUTS,
                             process_outputs=proc.PROCESS_OUTPUTS,
                             abstract=proc.ABSTRACT, hidden=proc.HIDDEN)

    def execute_workflow(self, identifier, context):
        process_task = base.get_process(identifier)

        start = tasks.job_started.s(context).set(**helpers.DEFAULT_QUEUE)

        workflow = process_task.s().set(**helpers.DEFAULT_QUEUE)

        succeeded = tasks.job_succeeded_workflow.s().set(**helpers.DEFAULT_QUEUE)

        canvas = start | workflow | succeeded;

        canvas.delay()

    def execute(self, identifier, variable, domain, operation, process, job, user):
        logger.info('Identifier %r', identifier)

        logger.info('Variable %r', variable)

        logger.info('Domain %r', domain)

        logger.info('Operation %r', operation)

        if identifier == 'CDAT.workflow':
            context = WorkflowOperationContext.from_data_inputs(variable,
                                                                domain,
                                                                operation)

            context.job = job

            context.user = user

            context.process = process

            self.execute_workflow(identifier, context)
        else:
            context = OperationContext.from_data_inputs(identifier, variable,
                                                        domain, operation)

            context.job = job

            context.user = user

            context.process = process

            process_func = base.get_process(identifier)

            if process_func.METADATA.get('inputs', 0) == 0:
                start = tasks.job_started.s(context).set(**helpers.DEFAULT_QUEUE)

                process = process_func.s().set(**helpers.DEFAULT_QUEUE)

                success = tasks.job_succeeded.s().set(**helpers.DEFAULT_QUEUE)

                canvas = start | process 

                canvas.delay()
            else:
                preprocess_chains = []

                for index in range(settings.WORKER_PER_USER):
                    filter = tasks.filter_inputs.s(index).set(**helpers.DEFAULT_QUEUE)

                    map = tasks.map_domain.s().set(**helpers.DEFAULT_QUEUE)

                    cache = tasks.check_cache.s().set(**helpers.DEFAULT_QUEUE)

                    chunks = tasks.generate_chunks.s().set(**helpers.DEFAULT_QUEUE)

                    preprocess_chains.append(celery.chain(filter, map, cache,
                                                          chunks))

                start = tasks.job_started.s(context).set(**helpers.DEFAULT_QUEUE)

                units = tasks.base_units.s().set(**helpers.DEFAULT_QUEUE)

                merge = tasks.merge.s().set(**helpers.DEFAULT_QUEUE)

                preprocess = start | units | celery.group(preprocess_chains) | merge

                process_chains = []

                for index in range(settings.WORKER_PER_USER):
                    ingress = tasks.ingress_chunk.s(index).set(**helpers.DEFAULT_QUEUE)

                    process = process_func.s(index).set(**helpers.DEFAULT_QUEUE)

                    process_chains.append(celery.chain(ingress, process))

                concat = tasks.concat.s().set(**helpers.DEFAULT_QUEUE)

                success = tasks.job_succeeded.s().set(**helpers.DEFAULT_QUEUE)

                cache = tasks.ingress_cache.s().set(**helpers.DEFAULT_QUEUE)

                cleanup = tasks.ingress_cleanup.s().set(**helpers.DEFAULT_QUEUE)

                finalize = concat | success | cache | cleanup

                canvas = preprocess | celery.group(process_chains) | finalize

                canvas.delay()
