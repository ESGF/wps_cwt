from __future__ import division
from future import standard_library
standard_library.install_aliases() # noqa
from builtins import str
from builtins import object
import os
import uuid
import coreapi
import requests

import cwt
from celery.utils.log import get_task_logger
from django.conf import settings

from wps import WPSError

logger = get_task_logger('wps.context')


class StateMixin(object):
    def __init__(self):
        self.job = None
        self.user = None
        self.process = None
        self.status = None
        self.client = None
        self.schema = None

    def init_state(self, data):
        self.job = data['job']

        self.user = data['user']

        self.process = data['process']

        self.init_api()

    def init_api(self):
        session = requests.Session()

        session.verify = not settings.DEBUG

        auth = coreapi.auth.BasicAuthentication(os.environ['API_USERNAME'], os.environ['API_PASSWORD'])

        transport = coreapi.transports.HTTPTransport(auth=auth, session=session)

        self.client = coreapi.Client(transports=[transport, ])

        schema_url = 'https://{!s}/internal_api/schema'.format(settings.INTERNAL_LB)

        self.schema = self.client.get(schema_url)

    def store_state(self):
        return {
            'job': self.job,
            'user': self.user,
            'process': self.process,
        }

    def action(self, keys, params=None):
        try:
            return self.client.action(self.schema, keys, params=params)
        except Exception as e:
            logger.exception(str(e))

            raise WPSError('Internal API call failed')

    def set_status(self, status, output=None, exception=None):
        params = {
            'job_pk': self.job,
            'status': status,
        }

        if output is not None:
            params['output'] = output

        if exception is not None:
            params['exception'] = exception

        output = self.action(['jobs', 'status', 'create'], params)

        self.status = output['id']

    def message(self, text, percent=None):
        if percent is None:
            percent = 0.0

        params = {
            'job_pk': self.job,
            'status_pk': self.status,
            'message': text,
            'percent': percent,
        }

        self.action(['jobs', 'status', 'message', 'create'], params)

    def accepted(self):
        self.set_status('ProcessAccepted')

    def started(self):
        self.set_status('ProcessStarted')

    def failed(self, exception):
        self.set_status('ProcessFailed', exception=exception)

    def succeeded(self, output):
        self.set_status('ProcessSucceeded', output=output)

    def processes(self):
        output = self.action(['process', 'list'])

        return output['results']

    def track_file(self, file):
        params = {
            'user_pk': self.user,
            'url': file.uri,
            'var_name': file.var_name,
        }

        self.action(['user', 'file', 'create'], params)


class WorkflowOperationContext(StateMixin, object):
    def __init__(self, inputs, domain, operation):
        self._inputs = inputs
        self.domain = domain
        self.operation = operation
        self.output = []

    @classmethod
    def from_data_inputs(cls, variable, domain, operation):
        for op in operation.values():
            inputs = []

            for x in op.inputs:
                if x in variable:
                    inputs.append(variable[x])
                elif x in operation:
                    inputs.append(operation[x])
                else:
                    raise WPSError('Unable to find input {!r}', x)

            op.inputs = inputs

            op.domain = domain.get(op.domain, None)

        instance = cls(variable, domain, operation)

        return instance

    @classmethod
    def from_dict(cls, data):
        variable = data['_inputs']

        domain = data['domain']

        operation = data['operation']

        for x in operation.values():
            inputs = []

            for y in x.inputs:
                if y in variable:
                    inputs.append(variable[y])
                elif y in operation:
                    inputs.append(operation[y])
                else:
                    raise WPSError('Failed to resolve input {!r}', y)

            x.inputs = inputs

            try:
                x.domain = domain[x.domain]
            except KeyError:
                pass

        instance = cls(variable, domain, operation)

        instance.output = [cwt.Variable.from_dict(x) for x in data['output']]

        instance.init_state(data)

        return instance

    def to_dict(self):
        data = {
            '_inputs': self._inputs,
            'domain': self.domain,
            'operation': self.operation,
            'output': self.output,
        }

        data.update(self.store_state())

        return data

    @property
    def inputs(self):
        return self._inputs.values()

    def output_ops(self):
        out_deg = dict((x, self.node_out_deg(x)) for x in self.operation.keys())

        return [self.operation[x] for x, y in out_deg.items() if y == 0]

    def node_in_deg(self, node):
        return len([x for x in node.inputs if x.name in self.operation])

    def node_out_deg(self, node):
        return len([x for x, y in self.operation.items() if any(node == z.name for z in y.inputs)])

    def find_neighbors(self, node):
        return [x for x, y in self.operation.items() if node in [x.name for x in y.inputs]]

    def topo_sort(self):
        in_deg = dict((x, self.node_in_deg(y)) for x, y in self.operation.items())

        neigh = dict((x, self.find_neighbors(x)) for x in in_deg.keys())

        queue = [x for x, y in in_deg.items() if y == 0]

        cnt = 0

        topo_order = []

        while queue:
            next = queue.pop(0)

            topo_order.append(self.operation[next])

            for x in neigh[next]:
                in_deg[x] -= 1

                if in_deg[x] == 0:
                    queue.append(x)

            cnt += 1

        if cnt != len(self.operation):
            raise WPSError('Failed to compute the graph')

        logger.info('Result of topo sort %r', topo_order)

        return topo_order

    def build_output_variable(self, var_name, name=None):
        local_path = self.generate_local_path()

        relpath = os.path.relpath(local_path, settings.WPS_PUBLIC_PATH)

        url = settings.WPS_DAP_URL.format(filename=relpath)

        self.output.append(cwt.Variable(url, var_name, name=name))

        return local_path

    def generate_local_path(self):
        filename = '{}.nc'.format(uuid.uuid4())

        base_path = os.path.join(settings.WPS_PUBLIC_PATH, str(self.user.id),
                                 str(self.job.id))

        if not os.path.exists(base_path):
            os.makedirs(base_path)

        return os.path.join(base_path, filename)


class OperationContext(StateMixin, object):
    def __init__(self, inputs=None, domain=None, operation=None):
        self.inputs = inputs
        self.domain = domain
        self.operation = operation
        self.output = []

    @classmethod
    def from_data_inputs(cls, identifier, variable, domain, operation):
        try:
            target_op = [x for x in list(operation.values()) if x.identifier ==
                         identifier][0]
        except IndexError:
            raise WPSError('Unable to find operation {!r}', identifier)

        logger.info('Target operation %r', target_op)

        target_domain = domain.get(target_op.domain, None)

        logger.info('Target domain %r', target_domain)

        try:
            target_inputs = [variable[x] for x in target_op.inputs]
        except KeyError as e:
            raise WPSError('Missing variable with name {!r}', str(e))

        logger.info('Target inputs %r', target_inputs)

        return cls(target_inputs, target_domain, target_op)

    @classmethod
    def from_dict(cls, data):
        try:
            inputs = [cwt.Variable.from_dict(x) for x in data['inputs']]

            if 'domain' in data and data['domain'] is not None:
                domain = cwt.Domain.from_dict(data['domain'])
            else:
                domain = None

            operation = cwt.Process.from_dict(data['operation'])
        except Exception:
            obj = cls()
        else:
            obj = cls(inputs, domain, operation)

        obj.output = [cwt.Variable.from_dict(x) for x in data['output']]

        obj.init_state(data)

        return obj

    def to_dict(self):
        data = {
            'inputs': self.inputs,
            'domain': self.domain,
            'operation': self.operation,
            'output': self.output,
        }

        data.update(self.store_state())

        return data

    @property
    def is_regrid(self):
        return 'gridder' in self.operation.parameters

    def build_output_variable(self, var_name, name=None):
        local_path = self.generate_local_path()

        relpath = os.path.relpath(local_path, settings.WPS_PUBLIC_PATH)

        url = settings.WPS_DAP_URL.format(filename=relpath)

        self.output.append(cwt.Variable(url, var_name, name=name))

        return local_path

    def generate_local_path(self):
        filename = '{}.nc'.format(uuid.uuid4())

        base_path = os.path.join(settings.WPS_PUBLIC_PATH, str(self.user.id),
                                 str(self.job.id))

        if not os.path.exists(base_path):
            os.makedirs(base_path)

        return os.path.join(base_path, filename)
