#! /usr/bin/env python

import types
from functools import partial

import cdms2
import cwt
import dask
import dask.array as da
from celery.utils.log import get_task_logger
from dask.distributed import Client
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from distributed.client import futures_of
from jinja2 import Environment, BaseLoader
from tornado.ioloop import IOLoop

from compute_tasks import base
from compute_tasks import context as ctx
from compute_tasks.context import state_mixin
from compute_tasks import managers
from compute_tasks import WPSError
from compute_tasks.dask_serialize import regrid_chunk

logger = get_task_logger('compute_tasks.cdat')


class DaskJobTracker(ProgressBar):
    def __init__(self, context, futures, scheduler=None, interval='100ms', complete=True):
        futures = futures_of(futures)

        context.message('Tracking {!r} futures', len(futures))

        super(DaskJobTracker, self).__init__(futures, scheduler, interval, complete)

        self.context = context
        self.last = None

        self.loop = IOLoop()

        loop_runner = LoopRunner(self.loop)
        loop_runner.run_sync(self.listen)

    def _draw_bar(self, **kwargs):
        logger.debug('_draw_bar %r', kwargs)

        remaining = kwargs.get('remaining', 0)

        all = kwargs.get('all', None)

        frac = (1 - remaining / all) if all else 1.0

        percent = int(100 * frac)

        logger.debug('Percent processed %r', percent)

        if self.last is None or self.last != percent:
            self.context.message('Processing', percent=percent)

            self.last = percent

    def _draw_stop(self, **kwargs):
        pass


def gather_workflow_outputs(context, interm, operations):
    delayed = []

    for output in operations:
        try:
            interm_value = interm.pop(output.name)
        except KeyError as e:
            raise WPSError('Failed to find intermediate {!s}', e)

        context.track_out_bytes(interm_value.nbytes)

        dataset = interm_value.to_xarray()

        output_name = '{!s}-{!s}'.format(output.name, output.identifier)

        local_path = context.build_output_variable(interm_value.var_name, name=output_name)

        context.message('Building output for {!r} - {!r}', output.name, output.identifier)

        logger.debug('Writing local output to %r', local_path)

        # Create an output file and store the future
        delayed.append(dataset.to_netcdf(local_path, compute=False))

    return delayed


def gather_inputs(identifier, fm, inputs):
    if identifier == 'CDAT.aggregate':
        gathered = [managers.InputManager.from_cwt_variables(fm, inputs)]
    else:
        gathered = [managers.InputManager.from_cwt_variable(fm, x) for x in inputs]

    return gathered


def build_workflow(fm, context):
    # Topologically sort the operations
    topo_order = context.topo_sort()

    # Hold the intermediate inputs
    interm = {}

    while topo_order:
        next = topo_order.pop(0)

        context.message('Processing operation {!r} - {!r}', next.name, next.identifier)

        if all(isinstance(x, cwt.Variable) for x in next.inputs):
            inputs = gather_inputs(next.identifier, fm, next.inputs)
        else:
            try:
                inputs = [interm[x.name] for x in next.inputs]
            except KeyError as e:
                raise WPSError('Missing intermediate data {!s}', e)

            # Create copy of single input, multiple inputs automatically create a new output
            if len(inputs) == 1:
                inputs = [inputs[0].copy(), ]

        process_func = PROCESS_FUNC_MAP[next.identifier]

        interm[next.name] = process_func(context, next, *inputs)

        context.message('Storing intermediate {!r}', next.name)

    return interm


WORKFLOW_ABSTRACT = """
This operation is used to store global values in workflows. Domain, regridders and parameters defined
here will become the default values on child operations.
"""


@base.register_process('CDAT', 'workflow', abstract=WORKFLOW_ABSTRACT, inputs='*')
@base.cwt_shared_task()
def workflow_func(self, context):
    """ Executes a workflow.

    Process a forest of operations. The graph can have multiple inputs and
    outputs.

    Args:
        context (WorkflowOperationContext): Current context.

    Returns:
        Updated context.
    """
    fm = managers.FileManager(context)

    interm = build_workflow(fm, context)

    context.message('Preparing to execute workflow')

    client = state_mixin.retry(8, 1)(Client)(context.extra['DASK_SCHEDULER'])

    try:
        delayed = []

        delayed.extend(gather_workflow_outputs(context, interm, context.output_ops()))

        delayed.extend(gather_workflow_outputs(context, interm, context.interm_ops()))

        context.message('Gathered {!s} outputs', len(delayed))

        with ctx.ProcessTimer(context):
            # Execute the futures
            fut = client.compute(delayed)

            # Track the progress
            DaskJobTracker(context, fut)
    except Exception as e:
        raise WPSError('Error executing process: {!r}', e)

    return context


def regrid(operation, *inputs):
    """ Build regridding from dask delayed.

    Converts the current dask array into an array of delayed functions. Next for each delayed function.
    another delayed function is linked. This function will regrid the chunk of data. Next the new delayed
    functions are converted back to dask arrays which are then concatenated back together. Lastly the
    spatial axes and bounds are updated.

    Args:
        operation: A cwt.Process instance.
        inputs: A list of compute_tasks.managers.FileManager instances.
    """
    gridder = operation.get_parameter('gridder', True)

    input = inputs[0]

    # Build a selector which is a dict mapping axis name to the desired range
    selector = dict((x, (input.axes[x][0], input.axes[x][-1])) for x in input.variable_axes)

    # Generate the target grid
    grid, tool, method = managers.generate_grid(gridder, selector), gridder.tool, gridder.method

    logger.info('Using grid shape %r tool %r method %r', grid.shape, tool, method)

    # Convert to delayed functions
    delayed = input.variable.to_delayed().squeeze()

    # Flatten the array
    if delayed.size == 1:
        logger.info('Flattening array')

        delayed = delayed.reshape((1, ))

    logger.info('Delayed shape %r', delayed.shape)

    # Create list of axis data that is present in the selector.
    axis_data = [input.axes[x] for x in selector.keys()]

    # Build list of delayed functions that will process the chunks
    regrid_delayed = [dask.delayed(regrid_chunk)(x, axis_data, grid, tool, method) for x in delayed]

    logger.info('Created regrid delayed functions')

    # Build a list of dask arrays created by the delayed functions
    regrid_arrays = [da.from_delayed(x, y.shape[:-2] + grid.shape, input.dtype)
                     for x, y in zip(regrid_delayed, input.blocks)]

    # Concatenated the arrays together
    input.variable = da.concatenate(regrid_arrays)

    logger.info('Concatenated %r arrays to output %r', len(regrid_arrays), input.variable)

    input.axes['lat'] = grid.getLatitude()

    input.vars['lat_bnds'] = cdms2.createVariable(
        input.axes['lat'].getBounds(),
        id='lat',
        axes=[input.axes[x] for x in input.vars_axes['lat_bnds']]
    )

    input.axes['lon'] = grid.getLongitude()

    input.vars['lon_bnds'] = cdms2.createVariable(
        input.axes['lon'].getBounds(),
        id='lon',
        axes=[input.axes[x] for x in input.vars_axes['lon_bnds']]
    )

    return input


FEAT_AXES = 'FEAT_AXES'
FEAT_CONST = 'FEAT_CONST'
FEAT_MULTI = 'FEAT_MULTI'


def process_input(context, operation, *inputs, **kwargs):
    """ Process inputs.

    Possible values for `supported`:
        single_input_axes: Process a single input array over some axes.
        single_input_constant: Process a single input array with a constant value.
        multiple_input: Process multiple inputs.

    Args:
        operation: A cwt.Process instance.
        inputs: A list of compute_tasks.managers.FileManager instances.
        process_func: An optional dask ufunc.
        supported: An optional dict defining the support processing types. See above for explanation.
    """
    for input in inputs:
        src, subset = input.subset(context.domain)

        context.track_src_bytes(src.nbytes)

        context.track_in_bytes(subset.nbytes)

        context.message('Data subset shape {!r}', input.shape)

    process_func = kwargs.get('process_func', None)

    if process_func is not None:
        axes = operation.get_parameter('axes')

        constant = operation.get_parameter('constant')

        logger.info('Axes %r Constant %r', axes, constant)

        if axes is not None:
            if not kwargs.get(FEAT_AXES, False):
                raise WPSError('Axes parameter is not supported by operation {!r}', operation.identifier)

            # Apply process to first input over axes
            output = process_single_input(axes.values, process_func, inputs[0])
        elif constant is not None:
            if not kwargs.get(FEAT_CONST, False):
                raise WPSError('Constant parameter is not supported by operation {!r}', operation.identifier)

            try:
                constant = float(constant.values[0])
            except ValueError:
                raise WPSError('Invalid constant value {!r} type {!s} expecting <class \'float\'>',
                               constant, type(constant))

            output = inputs[0]

            # Apply the process to the existing dask array
            output.variable = process_func(output.variable, constant)

            logger.info('Process output %r', output.variable)
        elif len(inputs) > 1:
            if not kwargs.get(FEAT_MULTI, False):
                raise WPSError('Multiple inputs are not supported by operation {!r}', operation.identifier)

            # Apply the process to all inputs
            output = process_multiple_input(process_func, *inputs)
        else:
            output = inputs[0]

            # Apply the process
            output.variable = process_func(output.variable)

            logger.info('Process output %r', output.variable)
    else:
        output = inputs[0]

    return output


def process_single_input(axes, process_func, input):
    """ Process single input.

    Args:
        axes: List of str axis names.
        process_func: A function with the same signature as process_input.
        input: A compute_tasks.managers.FileManager instance.

    Returns:
        A new instance of compute_tasks.managers.FileManager containing the results.
    """
    # Get present in the current data
    map_keys = list(input.variable_axes)

    # Reduce them from names to indices
    indices = tuple(map_keys.index(x) for x in axes)

    logger.info('Mapped axes %r to indices %r', axes, indices)

    # Apply the dask ufunc
    input.variable = process_func(input.variable, axis=indices)

    logger.info('Process output %r', input.variable)

    # Remove the axes that have been squashed
    for axis in axes:
        logger.info('Removing axis %r', axis)

        input.remove_axis(axis)

    return input


def process_multiple_input(process_func, input1, input2):
    """ Process multiple inputs.

    Note that some ufuncs will only process a single array over some dimensions,
    in this the function will need to be registed in the REQUIRES_STACK variable.

    Args:
        process_func: A function with the same signature as process_input.
        input1: A compute_tasks.managers.FileManager instance.
        input2: A compute_tasks.managers.FileManager instance.

    Returns:
        A new instance of compute_tasks.managers.FileManager containing the results.
    """
    # Create a copy to retain any metadata
    new_input = input1.copy()

    # Check if we need to stack the arrays.
    if process_func.__name__ in REQUIRES_STACK:
        stacked = da.stack([input1.variable, input2.variable])

        logger.info('Stacking inputs %r', stacked)

        new_input.variable = process_func(stacked, axis=0)

    else:
        new_input.variable = process_func(input1.variable, input2.variable)

    logger.info('Process output %r', new_input.variable)

    return new_input


# TODO remove this and utilize **kwargs on process_input to configure when stacking is required.
REQUIRES_STACK = [
    'sum_func',
    'max_func',
    'min_func',
    'average_func',
]

# Process descriptions used in abstracts.
DESCRIPTION_MAP = {
    'CDAT.abs': 'Computes the element-wise absolute value.',
    'CDAT.add': 'Adds an element-wise constant or another input.',
    'CDAT.aggregate': 'Aggregates multiple files over a temporal axis.',
    'CDAT.average': 'Computes the average over a set of axes or inputs.',
    'CDAT.divide': 'Divides element-wise by a constant or second input.',
    'CDAT.exp': 'Computes element-wise exponential.',
    'CDAT.log': 'Computes element-wise natural log',
    'CDAT.max': 'Computes the maximum over a set of axes, between a constant or a second input.',
    'CDAT.min': 'Computes the minimum over a set of axes, between a constant or a second input.',
    'CDAT.multiply': 'Multiplies element-wise by a constant or a second input.',
    'CDAT.power': 'Computes the element-wise power by a constant.',
    'CDAT.regrid': 'Regrids input by target grid.',
    'CDAT.subset': 'Subsets an input to desired domain.',
    'CDAT.subtract': 'Subtracts element-wise constant or a second input.',
    'CDAT.sum': 'Computes the sum over a set of axes.',
}

"""
This dict maps identifiers to functions that will be used for processing. The function ``process_input`` is
the main entrypoint. Partials are used here to provide some configuration for how each process is executed.
See the doctstring of ``process_input`` for additional details.
"""
PROCESS_FUNC_MAP = {
    'CDAT.abs': partial(process_input, process_func=da.absolute),
    'CDAT.add': partial(process_input, process_func=da.add, FEAT_CONST=True, FEAT_MULTI=True),
    'CDAT.aggregate': process_input,
    'CDAT.average': partial(process_input, process_func=da.average, FEAT_MULTI=True),
    'CDAT.divide': partial(process_input, process_func=da.divide, FEAT_CONST=True, FEAT_MULTI=True),
    # Disabled due to overflow issue, setting dtype=float64 works for dask portion but xarray is writing Inf.
    # 'CDAT.exp': partial(process_input, process_func=da.exp),
    'CDAT.log': partial(process_input, process_func=da.log),
    'CDAT.max': partial(process_input, process_func=da.max, FEAT_AXES=True, FEAT_CONST=True, FEAT_MULTI=True),
    'CDAT.min': partial(process_input, process_func=da.min, FEAT_AXES=True, FEAT_CONST=True, FEAT_MULTI=True),
    'CDAT.multiply': partial(process_input, process_func=da.multiply, FEAT_CONST=True, FEAT_MULTI=True),
    'CDAT.power': partial(process_input, process_func=da.power, FEAT_CONST=True),
    # 'CDAT.regrid': regrid,
    'CDAT.subset': process_input,
    'CDAT.subtract': partial(process_input, process_func=da.subtract, FEAT_CONST=True, FEAT_MULTI=True),
    'CDAT.sum': partial(process_input, process_func=da.sum, FEAT_AXES=True),
}


def render_abstract(description, func, template):
    """ Renders an abstract for a process.

    This function will use a jinja2 template and render out an abstract for a process. The keyword arguments
    for the ``func`` function are used to enable details in the abstract.

    Args:
        description (str): The process description.
        func (function): The process function.
        template (jinja2.Template): The jinja2 template that will be used to render the abstract.

    Returns:
        str: The abstract as a string.
    """
    axes = False
    const = False
    multi = False

    try:
        kwargs = func.keywords
    except AttributeError:
        pass
    else:
        axes = FEAT_AXES in kwargs and kwargs[FEAT_AXES]

        const = FEAT_CONST in kwargs and kwargs[FEAT_CONST]

        multi = FEAT_MULTI in kwargs and kwargs[FEAT_MULTI]

    return template.render(description=description, axes=axes, const=const, multi=multi)


BASE_ABSTRACT = """
{{ description }}
{%- if multi %}

Supports multiple inputs.
{%- endif %}
{%- if (axes or constant) %}

Optional parameters:
{%- if axes %}
    axes: A list of axes to operate on. Multiple values should be separated by "|" e.g. "lat|lon".
{%- endif %}
{%- if constant %}
    constant: An integer or float value that will be applied element-wise.
{%- endif %}
{%- endif %}
"""


def process_wrapper(self, context):
    """ Wrapper function for a process.

    This function acts as the main entrypoint of a Celery task. It represents a single process e.g. Subset, Aggregate,
    etc. The function calls ``workflow_func`` since a single process is just a workflow with a single task.

    Args:
        context (OperationContext): The OperationContext holding all details of the current job.
    """
    return workflow_func(context)


def copy_function(f, operation):
    """ Creates a unique version of a function.

    Copies function ``f`` giving it a unique name using ``operation``.

    Args:
        f (FunctionType): The function to copy.
        operation (str): The unique identifier for the function.

    Returns:
        FunctionType: A new function.
    """
    name = '{!s}_func'.format(operation)

    return types.FunctionType(f.__code__, f.__globals__, name=name, argdefs=f.__defaults__, closure=f.__closure__)


def discover_processes():
    """ Discovers and binds functions to `cdat` module.

    This function iterates over PROCESS_FUNC_MAP, generating a description, creating a Celery task, registering it with
    the backend and binding it to the "cdat" module.

    Returns:
        list: List of dict, describing each registered process.
    """
    from compute_tasks import base
    from compute_tasks import cdat

    # Use jinja2 to template process abstract
    template = Environment(loader=BaseLoader).from_string(BASE_ABSTRACT)

    for name, func in PROCESS_FUNC_MAP.items():
        module, operation = name.split('.')

        # Create a unique function for each process
        p = copy_function(process_wrapper, operation)

        # Decorate the new function as a Celery task
        shared = base.cwt_shared_task()(p)

        # Render the abstract
        abstract = render_abstract(DESCRIPTION_MAP[name], func, template)

        inputs = 1

        try:
            if FEAT_MULTI in func.keywords:
                inputs = 2
            elif name == 'CDAT.aggregate':
                inputs = '*'
        except AttributeError:
            # Handle where no func is supplied
            pass

        # Decorate the Celery task as a registered process
        register = base.register_process(module, operation, abstract=abstract, inputs=inputs)(shared)

        # Bind the new function to the "cdat" module
        setattr(cdat, p.__name__, register)

    return base.REGISTRY.values()
