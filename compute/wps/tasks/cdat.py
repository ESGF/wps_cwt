#! /usr/bin/env python

import datetime
import json
import os
import re
import uuid

import cdms2
import cwt
import cdutil
from cdms2 import MV2 as MV
from celery.task.control import inspect
from celery.utils.log import get_task_logger
from django.conf import settings
from django.utils import timezone

from wps import helpers
from wps import metrics
from wps import models
from wps import WPSError
from wps.tasks import base

logger = get_task_logger('wps.tasks.cdat')

OUTPUT = cwt.wps.process_output_description('output', 'output', 'application/json')

PATTERN_AXES_REQ = 'CDAT\.(min|max|average|sum)'

def retrieve_data(infile, outfile, var_name, grid, gridder, base_units, mapped=None):
    """ Retrieves data and writes to output.

    Reads a subset of the input data where mapped is a dict selector. Will
    regrid the data and rebase the time access. Finally write data to output
    file.

    Args:
        infile: A cdms2.dataset.CdmsFile object.
        outfile: A cdms2.dataset.CdmsFile object.
        var_name: A str variable name.
        grid: A cdms2.grid.* object.
        gridder: A cwt.Gridder object.
        base_units: A str units to rebase the time axis.
        mapped: A dict selector to subset data.

    Returns:
        None

    """
    if mapped is None:
        mapped = {}

    logger.info('Retrieving data selector %r, grid %r, base_units %r', mapped,
                grid, base_units)

    data = infile(var_name, **mapped)

    if grid is not None:
        data = data.regrid(grid, regridTool=gridder.tool, regridMethod=gridder.method)

    if base_units is not None:
        data.getTime().toRelativeTime(str(base_units))

    outfile.write(data, id=var_name)

    return data

def retrieve_data_cached(self, infile, outfile, var_name, grid, gridder, base_units, mapped, chunk_axis, chunk_list, **kwargs):
    """ Retrieves cached data.

    Really just a convenience method. Splits out mapped, chunk_axis and 
    chunk_list from a dict. 

    Will subset a grid if supplied and then retrieves each chunk from the cached
    file.

    """
    if grid is not None:
        grid = self.subset_grid(grid, mapped)

    for chunk in chunk_list:
        mapped.update({chunk_axis: chunk})

        start = self.get_now()

        data = retrieve_data(infile, outfile, var_name, grid, gridder, base_units, mapped)

        elapsed = self.get_now() - start

        metrics.CACHE_BYTES.inc(data.nbytes)

        metrics.CACHE_SECONDS.inc(elapsed.total_seconds())

def base_retrieve(self, attrs, keys, operation, var_name, base_units, output_path, job_id):
    """ Retrieve file(s).

    This is the base for aggregate, subset and regrid.

    Expected format for "attrs" argument.

    {
        "key": {
            "path": "https://aims3.llnl.gov/path/filename.nc",
            "chunked_axis": "time",
            "chunks": {
                "time": [slice(0, 10), slice(10, 12)],
            },
            "mapped": {
                "time": slice(0, 10),
                "lat": slice(0, 100),
                "lon": slice(0, 200),
            },
        }
    }

    Args:
        attrs: A list of dict or dict from previous tasks.
        cached: A list of dict of cached portions.
        operation: A cwt.Process object.
        var_name: A str variable name.
        base_units: A str base_units to be used.
        job_id: An int of the current job id.

    Returns:
        The input attrs.
    """
    job = self.load_job(job_id)

    gridder = operation.get_parameter('gridder')

    if not isinstance(attrs, list):
        attrs = [attrs]

    attrs = dict(y for x in attrs for y in x.items())

    grid = None
    selector = None

    start = self.get_now()

    try:
        os.makedirs(os.path.dirname(output_path))
    except OSError:
        raise WPSError('Failed to create output directory')

    with self.open(output_path, 'w') as outfile:
        # Expect the keys to be given in a sortable format
        for key in sorted(keys):
            current = attrs[key]

            with self.open(current['path']) as infile:
                # Generate the grid once
                if grid is None and gridder is not None:
                    grid = self.generate_grid(gridder)

                    self.update(job, 'Generated grid {!r}', grid)

                if 'cached' in current:
                    retrieve_data_cached(self, infile, outfile, var_name, grid, gridder, base_units, **current)

                    self.update(job, 'Building file from cache')
                else:
                    # Subset the grid to the target shape
                    if selector is None and grid is not None:
                        selector = self.generate_selector(infile[var_name])

                        self.update(job, 'Generated subset selector {!r}',
                                    selector)

                        grid = self.subset_grid(grid, selector)

                        self.update(job, 'Subsetting grid {!r}', grid)

                    retrieve_data(infile, outfile, var_name, grid, gridder, base_units)

                    self.update(job, 'Building file from ingressed data {}',
                                current['path'].split('/')[-1])

    self.update(job, 'Finished building file {}',
                output_path.split('/')[-1])

    elapsed = self.get_now() - start

    stat = os.stat(output_path)

    metrics.PROCESS_BYTES.labels(operation.identifier).inc(stat.st_size)

    metrics.PROCESS_SECONDS.labels(operation.identifier).inc(elapsed.total_seconds())

    return attrs

def base_process(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id):
    """ Process file.

    Expected format for attrs argument.

    {
        "key": {
            "path": "file:///path/filename.nc",
            "mapped": {
                "time": slice(0, 200, 1),
                ...
            }
        }
    }

    Args:
        attrs: A dict describing the file to be processed.
        key: A key used to identify the data in attrs.
        operation: A cwt.Process object.
        var_name: A str variable name.
        base_units: A str units used to base the time axis.
        axes: A list of str names of axes to operate over.
        output_path: A str containing the output path.
        job_id: An int of the current job id. 

    Returns:
        The input attrs dict.
    """
    job = self.load_job(job_id)

    self.update(job, 'Processing chunk of "{}" with {} over {}', var_name,
                operation.identifier, axes)

    inp = attrs[key]

    gridder = operation.get_parameter('gridder')

    weightoptions = operation.get_parameter('weightoptions')

    mapped = inp.get('mapped', {})

    # Read input data
    with self.open(inp['path']) as infile:
        start = self.get_now()

        data = infile(var_name, **mapped)

        elapsed = self.get_now() - start

        metrics.CACHE_BYTES.inc(data.nbytes)

        metrics.CACHE_SECONDS.inc(elapsed.total_seconds())

    start = self.get_now()

    # Generate grid if needed
    if gridder is not None:
        # If we're processing an ingressed file we need to generate a selector
        # to subset the target grid.
        if len(mapped) == 0:
            mapped = self.generate_selector(data)

        grid = self.generate_grid(gridder)

        grid = self.subset_grid(grid, mapped)

        self.update(job, 'Regridding to {}', grid.shape)

        data = data.regrid(grid, regridTool=gridder.tool, regridMethod=gridder.method)

    # Grab the indexes of the axes
    axes_index = [data.getAxisIndex(str(x)) for x in axes]

    if operation.identifier == 'CDAT.average':
        axes_sig = ''.join([str(x) for x in axes_index])

        if weightoptions is not None:
            weightoptions = weightoptions.values[0]

        logger.info('Averaging over axes %r', axes_sig)

        data = self.PROCESS(data, axis=axes_sig, weights=weightoptions)
    else:
        # Process over all axes except the chunking axis
        for axis in axes_index:
            logger.info('Processing %r over axis %r', self.PROCESS, axis)

            data = self.PROCESS(data, axis=axis)

    with self.open(output_path, 'w') as outfile:
        outfile.write(data, id=var_name)

    elapsed = self.get_now() - start

    self.update(job, 'Finished processing chunk')

    attrs[key]['elapsed'] = elapsed.total_seconds()

    return attrs

@base.cwt_shared_task()
def concat_process_output(self, attrs, input_paths, operation, var_name, chunked_axis, output_path, job_id):
    """ Concatenates inputs over chunked_axis.

    Args:
        attrs: A dict or list of dicts from previous tasks.
        input_paths: A list of inputs to be concatenated.
        var_name: A str variable name.
        chunked_axis: A str name of the chunked axis.
        output_path: A str path to write the output to.
        job_id: An int referencing the associated job.

    Returns:
        The input attrs value.
    """
    job = self.load_job(job_id)

    self.update(job, 'Concatenating {} chunks over {}', len(input_paths),
                chunked_axis)

    data_list = []

    start = self.get_now()

    for input_path in sorted(input_paths):
        with self.open(input_path) as infile:
            data_list.append(infile(var_name))

    axis_index = data_list[0].getAxisIndex(chunked_axis)

    data = MV.concatenate(data_list, axis=axis_index)

    with self.open(output_path, 'w') as outfile:
        outfile.write(data, id=var_name)

    elapsed = self.get_now() - start

    new_attrs = {}

    for item in attrs:
        new_attrs.update(item)

    self.update(job, 'Finished concatentating chunks, final shape {}',
                data.shape)

    stat = os.stat(output_path)

    metrics.PROCESS_BYTES.labels(operation.identifier).inc(stat.st_size)

    total_seconds = sum(x['elapsed'] for x in new_attrs.values())

    total_seconds += elapsed.total_seconds()

    metrics.PROCESS_SECONDS.labels(operation.identifier).inc(total_seconds)

    return new_attrs

SNG_DATASET_SNG_INPUT = {
    'datasets': 1,
    'inputs': 1,
}

SNG_DATASET_MULTI_INPUT = {
    'datasets': 1,
    'inputs': '*',
}

@base.register_process('CDAT.regrid', abstract="""
                       Regrids a variable to designated grid. Required parameter named "gridder".
                       """, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def regrid(self, attrs, keys, operation, var_name, base_units, output_path, job_id):
    return base_retrieve(self, attrs, keys, operation, var_name, base_units, output_path, job_id)

@base.register_process('CDAT.subset', abstract="""
                       Subset a variable by provided domain. Supports regridding.
                       """, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def subset(self, attrs, keys, operation, var_name, base_units, output_path, job_id):
    return base_retrieve(self, attrs, keys, operation, var_name, base_units, output_path, job_id)

@base.register_process('CDAT.aggregate', abstract="""
                       Aggregate a variable over multiple files. Supports subsetting and regridding.
                       """, metadata=SNG_DATASET_MULTI_INPUT)
@base.cwt_shared_task()
def aggregate(self, attrs, keys, operation, var_name, base_units, output_path, job_id):
    return base_retrieve(self, attrs, keys, operation, var_name, base_units, output_path, job_id)

@base.register_process('CDAT.average', abstract=""" 
                       Computes the average over axes. 

                       Required parameters:
                         axes: A list of axes to operate on. Should be separated by "|".

                       Optional parameters:
                         weightoptions: A string whos value is "generate",
                           "equal", "weighted", "unweighted". See documentation
                           at https://cdat.llnl.gov/documentation/utilities/utilities-1.html
                       """, process=cdutil.averager, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def average(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id):
    return base_process(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id)

@base.register_process('CDAT.sum', abstract=""" 
                       Computes the sum over an axis. Requires singular parameter named "axes" 
                       whose value will be used to process over. The value should be a "|" delimited
                       string e.g. 'lat|lon'.
                       """, process=MV.sum, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def summation(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id):
    return base_process(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id)

@base.register_process('CDAT.max', abstract=""" 
                       Computes the maximum over an axis. Requires singular parameter named "axes" 
                       whose value will be used to process over. The value should be a "|" delimited
                       string e.g. 'lat|lon'.
                       """, process=MV.max, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def maximum(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id):
    return base_process(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id)

@base.register_process('CDAT.min', abstract="""
                       Computes the minimum over an axis. Requires singular parameter named "axes" 
                       whose value will be used to process over. The value should be a "|" delimited
                       string e.g. 'lat|lon'.
                       """, process=MV.min, metadata=SNG_DATASET_SNG_INPUT)
@base.cwt_shared_task()
def minimum(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id):
    return base_process(self, attrs, key, operation, var_name, base_units, axes, output_path, job_id)
