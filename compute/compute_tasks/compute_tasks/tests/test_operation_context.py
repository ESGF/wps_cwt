import json

import cwt
import pytest

from compute_tasks import WPSError
from compute_tasks.context import operation


class CWTData(object):
    def __init__(self):
        self.v1 = cwt.Variable('file:///test1.nc', 'tas')
        self.v2 = cwt.Variable('file:///test2.nc', 'tas')
        self.v3 = cwt.Variable('file:///test3.nc', 'pr')
        self.v4 = cwt.Variable('file:///test3.nc', 'prw')

        self.d0 = cwt.Domain(time=slice(10, 20, 2))

        self.op1 = cwt.Process(identifier='CDAT.aggregate')
        self.op2 = cwt.Process(identifier='CDAT.workflow')

        self.aggregate = self.op1
        self.workflow = self.op2
        self.max = cwt.Process(identifier='CDAT.max')
        self.min = cwt.Process(identifier='CDAT.min')
        self.subtract = cwt.Process(identifier='CDAT.subtract')
        self.divide = cwt.Process(identifier='CDAT.divide')
        self.std = cwt.Process(identifier='CDAT.std')

    def sample_workflow(self):
        self.std.add_inputs(self.aggregate)
        self.subtract.add_inputs(self.aggregate)
        self.divide.add_inputs(self.subtract, self.std)
        self.max.add_inputs(self.aggregate)
        self.min.add_inputs(self.aggregate)
        self.workflow.add_inputs(self.divide, self.max, self.min)

        data_inputs = {
            'variable': [],
            'domain': [],
            'operation': [
                self.aggregate.to_dict(),
                self.std.to_dict(),
                self.subtract.to_dict(),
                self.divide.to_dict(),
                self.max.to_dict(),
                self.min.to_dict(),
                self.workflow.to_dict(),
            ],
        }

        return data_inputs

    def process_mixed_inputs(self):
        s1 = cwt.Process(identifier='CDAT.subset', name='s1')
        s1.add_inputs(self.v3)

        m = cwt.Process(identifier='CDAT.merge', name='m')
        m.add_inputs(s1, self.v4)

        return {
            'variable': [self.v3.to_dict(), self.v4.to_dict()],
            'domain': [],
            'operation': [s1.to_dict(), m.to_dict()],
        }

    def workflow_with_rename(self):
        s1 = cwt.Process(identifier='CDAT.subset', name='s1')
        s1.add_inputs(self.v3)

        s2 = cwt.Process(identifier='CDAT.subset', name='s2')
        s2.add_inputs(self.v4)

        m = cwt.Process(identifier='CDAT.merge', name='m')
        m.add_inputs(s1, s2)

        c = cwt.Process(identifier='CDAT.count', name='c')
        c.add_inputs(m)
        c.add_parameters(rename=['pr', 'pr_count'])

        s = cwt.Process(identifier='CDAT.sum', name='s')
        s.add_inputs(c)
        s.add_parameters(variable=['pr_count'], rename=['pr_count', 'pr_sum'])

        w = cwt.Process(identifier='CDAT.workflow')
        w.add_inputs(s)

        return {
            'variable': [self.v3.to_dict(), self.v4.to_dict()],
            'domain': [],
            'operation': [s1.to_dict(), s2.to_dict(), m.to_dict(), c.to_dict(), s.to_dict(), w.to_dict()]
        }

@pytest.fixture(scope='function')
def cwt_data():
    return CWTData()


def test_topo_sort_mixed_inputs(cwt_data):
    workflow = cwt_data.process_mixed_inputs()

    ctx = operation.OperationContext.from_data_inputs('CDAT.merge', workflow)

    output = dict((x.name, y) for x, y in ctx.topo_sort())


def test_topo_sort(cwt_data):
    workflow = cwt_data.workflow_with_rename()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    output = dict((x.name, y) for x, y in ctx.topo_sort())

    assert 's1' in output
    assert output['s1'] == ['pr']
    assert 's2' in output
    assert output['s2'] == ['prw']
    assert 'm' in output
    assert sorted(output['m']) == ['pr', 'prw']
    assert 'c' in output
    assert sorted(output['c']) == ['pr', 'prw']
    assert 's' in output
    assert sorted(output['s']) == ['pr_count', 'prw']


def test_find_neighbors(cwt_data):
    workflow = cwt_data.sample_workflow()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    ops = ctx.interm_ops()

    node = [x for x in ops if x.identifier == 'CDAT.subtract'][0]

    neigh = ctx.find_neighbors(node.name)

    assert len(neigh) == 1
    assert neigh[0] == cwt_data.divide.name


def test_node_out_deg(cwt_data):
    workflow = cwt_data.sample_workflow()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    ops = ctx.interm_ops()

    node = [x for x in ops if x.identifier == 'CDAT.subtract'][0]

    assert ctx.node_out_deg(node) == 1


def test_node_in_deg(cwt_data):
    workflow = cwt_data.sample_workflow()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    ops = ctx.interm_ops()

    node = [x for x in ops if x.identifier == 'CDAT.subtract'][0]

    assert ctx.node_in_deg(node) == 1


def test_interm_ops(cwt_data):
    workflow = cwt_data.sample_workflow()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    ops = ctx.interm_ops()

    assert len(ops) == 3
    assert set([x.identifier for x in ops]) == set(['CDAT.aggregate', 'CDAT.std', 'CDAT.subtract'])


def test_output_ops(cwt_data):
    workflow = cwt_data.sample_workflow()

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', workflow)

    ops = ctx.output_ops()

    assert len(ops) == 3
    assert set([x.identifier for x in ops]) == set(['CDAT.divide', 'CDAT.min', 'CDAT.max'])


def test_to_dict(cwt_data):
    data = {
        '_variable': {
            'v1': cwt_data.v1,
        },
        '_domain': {
            'd0': cwt_data.d0,
        },
        '_operation': {
            'op1': cwt_data.op1,
        },
        'output': [
            cwt_data.v1,
        ],
        'extra': {},
        'job': 0,
        'user': 0,
        'process': 0,
        'gdomain': None,
        'gparameters': {},
        '_sorted': [],
        'input_var_names': {},
    }

    ctx = operation.OperationContext.from_dict(data)

    output = ctx.to_dict()

    assert '_variable' in output
    assert '_domain' in output
    assert '_operation' in output
    assert 'gdomain' in output
    assert 'gparameters' in output
    assert 'output' in output


def test_from_dict(cwt_data):
    data = {
        '_variable': {
            'v1': cwt_data.v1,
        },
        '_domain': {
            'd0': cwt_data.d0,
        },
        '_operation': {
            'op1': cwt_data.op1,
        },
        'output': [
            cwt_data.v1,
        ],
        'extra': {},
        'job': 0,
        'user': 0,
        'process': 0,
        'operation': cwt_data.op1,
        'gdomain': None,
        'gparameters': {},
        '_sorted': [],
        'input_var_names': {},
    }

    ctx = operation.OperationContext.from_dict(data)

    assert len(ctx._variable) == 1
    assert len(ctx._domain) == 1
    assert len(ctx._operation) == 1
    assert len(ctx.output) == 1


def test_from_data_inputs_missing_operation():
    data_inputs = {
        'variable': [],
        'domain': [],
        'operation': [],
    }

    with pytest.raises(WPSError):
        operation.OperationContext.from_data_inputs('CDAT.subset', data_inputs)


def test_from_data_inputs_invalid_input_workflow(cwt_data):
    cwt_data.op1.add_inputs(cwt_data.v1, cwt_data.v2)

    cwt_data.op2.add_inputs(cwt_data.op1)

    data_inputs = {
        'variable': [],
        'domain': [],
        'operation': [
            cwt_data.op1.to_dict(),
            cwt_data.op2.to_dict(),
        ],
    }

    with pytest.raises(WPSError):
        operation.OperationContext.from_data_inputs('CDAT.workflow', data_inputs)


def test_from_data_inputs_global(cwt_data):
    cwt_data.op1.add_inputs(cwt_data.v1, cwt_data.v2)
    cwt_data.op1.add_parameters(constant='10')

    cwt_data.op2.add_inputs(cwt_data.op1)
    cwt_data.op2.domain = cwt_data.d0
    cwt_data.op2.add_parameters(axes=['time'])

    data_inputs = {
        'variable': [
            cwt_data.v1.to_dict(),
            cwt_data.v2.to_dict(),
        ],
        'domain': [
            cwt_data.d0.to_dict(),
        ],
        'operation': [
            cwt_data.op1.to_dict(),
            cwt_data.op2.to_dict(),
        ],
    }

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', data_inputs)

    assert len(ctx._variable) == 2
    assert len(ctx._domain) == 1
    assert len(ctx._operation) == 1

    assert ctx.gdomain.name == cwt_data.d0.name
    assert ctx.gparameters.keys() == cwt_data.op2.parameters.keys()


def test_from_data_inputs_workflow(cwt_data):
    cwt_data.op1.add_inputs(cwt_data.v1, cwt_data.v2)
    cwt_data.op1.domain = cwt_data.d0

    cwt_data.op2.add_inputs(cwt_data.op1)

    data_inputs = {
        'variable': [
            cwt_data.v1.to_dict(),
            cwt_data.v2.to_dict(),
        ],
        'domain': [
            cwt_data.d0.to_dict(),
        ],
        'operation': [
            cwt_data.op1.to_dict(),
            cwt_data.op2.to_dict(),
        ],
    }

    ctx = operation.OperationContext.from_data_inputs('CDAT.workflow', data_inputs)

    assert len(ctx._variable) == 2
    assert len(ctx._domain) == 1
    assert len(ctx._operation) == 1


def test_from_data_inputs(cwt_data):
    cwt_data.op1.add_inputs(cwt_data.v1, cwt_data.v2)
    cwt_data.op1.domain = cwt_data.d0

    data_inputs = {
        'variable': [
            cwt_data.v1.to_dict(),
            cwt_data.v2.to_dict(),
        ],
        'domain': [
            cwt_data.d0.to_dict(),
        ],
        'operation': [
            cwt_data.op1.to_dict(),
        ],
    }

    ctx = operation.OperationContext.from_data_inputs('CDAT.aggregate', data_inputs)

    assert len(ctx._variable) == 2
    assert len(ctx._domain) == 1
    assert len(ctx._operation) == 1


def test_resolve_dependencies_globals(cwt_data):
    variable = {
        'v1': cwt_data.v1,
    }

    domain = {}

    cwt_data.op1.add_inputs('v1')

    operation_ = {
        'op1': cwt_data.op1,
    }

    parameters = {
        'axes': cwt.NamedParameter('axes', 'time'),
    }

    operation.OperationContext.resolve_dependencies(variable, domain, operation_, cwt_data.d0, parameters)

    assert cwt_data.op1.inputs[0] == cwt_data.v1
    assert cwt_data.op1.domain == cwt_data.d0
    assert len(cwt_data.op1.parameters) == 1


def test_resolve_dependencies_input_error(cwt_data):
    variable = {
        'v1': cwt_data.v1,
    }

    domain = {
        'd0': cwt_data.d0,
    }

    cwt_data.op1.add_inputs('v1', 'v2')
    cwt_data.op1.domain = 'd0'

    operation_ = {
        'op1': cwt_data.op1,
    }

    with pytest.raises(WPSError):
        operation.OperationContext.resolve_dependencies(variable, domain, operation_)


def test_resolve_dependencies(cwt_data):
    variable = {
        'v1': cwt_data.v1,
    }

    domain = {
        'd0': cwt_data.d0,
    }

    cwt_data.op1.add_inputs('v1')
    cwt_data.op1.domain = 'd0'

    operation_ = {
        'op1': cwt_data.op1,
    }

    operation.OperationContext.resolve_dependencies(variable, domain, operation_)

    assert cwt_data.op1.inputs[0] == cwt_data.v1
    assert cwt_data.op1.domain == cwt_data.d0


def test_decode_data_inputs_json_exception(cwt_data):
    data_inputs = {
        'domain': json.dumps([
            cwt_data.d0.to_dict(),
        ]),
        'operation': json.dumps([
            cwt_data.op1.to_dict(),
        ]),
    }

    with pytest.raises(WPSError):
        operation.OperationContext.decode_data_inputs(data_inputs)


def test_decode_data_inputs(cwt_data):
    data_inputs = {
        'variable': [
            cwt_data.v1.to_dict(),
            cwt_data.v2.to_dict(),
        ],
        'domain': [
            cwt_data.d0.to_dict(),
        ],
        'operation': [
            cwt_data.op1.to_dict(),
        ],
    }

    output = operation.OperationContext.decode_data_inputs(data_inputs)

    assert len(output[0]) == 2
    assert len(output[1]) == 1
    assert len(output[2]) == 1
