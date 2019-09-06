import re
import logging
from collections import OrderedDict

import cdms2
import cwt
import dask
import dask.array as da
import xarray as xr

from compute_tasks import metrics_ as metrics
from compute_tasks import WPSError
from compute_tasks.dask_serialize import retrieve_chunk

logger = logging.getLogger('compute_tasks.manager')


def format_dimension(dim):
    """ Formats a dimension.

    Formats a cwt.Dimension to either a slice or tuple. A slice denotes that the axis
    will not need to be mapped, whereas a tuple will need mapping.

    Args:
        dim: A cwt.Dimension to format.

    Returns:
        A slice or tuple.
    """

    if dim.crs == cwt.VALUES:
        data = (dim.start, dim.end, dim.step)
    elif dim.crs == cwt.INDICES:
        data = slice(dim.start, dim.end, dim.step)
    elif dim.crs == cwt.TIMESTAMPS:
        data = (dim.start, dim.end, dim.step)
    else:
        raise WPSError('Unknown dimension CRS %r', dim)

    logger.info('Formatted dimension as %r', data)

    return data


def domain_to_dict(domain):
    """ Converts a domain to a dict.

    Args:
        domain: A cwt.Domain to convert.

    Returns:
        A dict mapping dimension names to the formatted dimension value.
    """
    if domain is None:
        return {}

    output = dict((x, format_dimension(y)) for x, y in domain.dimensions.items())

    logger.info('Converted domain to %r', output)

    return output


def map_dimension(dim, axis):
    if not isinstance(dim, slice):
        try:
            interval = axis.mapInterval(dim[:2])
        except TypeError:
            mapped = None
        else:
            if len(dim) > 2:
                step = dim[2]
            else:
                step = 1

            # Prevent wrapping by taking min of stop value and len of axis
            mapped = slice(interval[0], min(interval[1], len(axis)), step)
    else:
        mapped = dim

    logger.info('Mapped dimension %r shape %r -> %r', axis.id, dim, mapped)

    return mapped


def map_domain(domain, axes):
    domain_dict = domain_to_dict(domain)

    logger.info('Mapping domain %r', domain_dict)

    map = OrderedDict()

    for name, value in axes.items():
        try:
            dim = domain_dict[name]
        except KeyError:
            map[name] = slice(None, None, None)
        else:
            map[name] = map_dimension(dim, value)

    return map


def subset_grid(grid, selector):
    target = cdms2.MV2.ones(grid.shape)

    logger.debug('Target grid %r', target.shape)

    target.setAxisList(grid.getAxisList())

    logger.info('Subsetting grid with selector %r', selector)

    target = target(**selector)

    logger.debug('Target grid new shape %r', target.shape)

    return target.getGrid()


def parse_uniform_arg(value):
    result = re.match(r'^(-?\d*\.?\d+):(\d*\.?\d+):(\d*\.?\d+)$', value)

    if not result:
        raise WPSError('Uniform grid argument {!r} does not match expected format '
                       '<start>:<number of gridlines>:<delta>', value)

    groups = result.groups()

    return float(groups[0]), float(groups[1]), float(groups[2])


def generate_user_defined_grid(gridder):
    try:
        grid_type, grid_param = gridder.grid.split('~')
    except AttributeError:
        return None
    except ValueError:
        raise WPSError('Error parsing grid definition {!r} expecting format <type>~<params>', gridder.grid)

    logger.info('Parsed grid type %r and parameters %r', grid_type, grid_param)

    if grid_type.lower() == 'uniform':
        result = re.match('^([^x]+)x(.+)$', grid_param)

        if result is None:
            raise WPSError('Error parsing uniform grid format {!r}', grid_param)

        start_lat, nlat, delta_lat = parse_uniform_arg(result.group(1))

        logger.info('Parsed Latitude start: %r nlat: %r delat: %r', start_lat, nlat, delta_lat)

        start_lon, nlon, delta_lon = parse_uniform_arg(result.group(2))

        logger.info('Parsed Longitude start: %r nlat: %r delat: %r', start_lon, nlon, delta_lon)

        grid = cdms2.createUniformGrid(start_lat, nlat, delta_lat, start_lon, nlon, delta_lon)

        logger.info('Created uniform grid {!r}', grid.shape)
    elif grid_type.lower() == 'gaussian':
        try:
            nlats = int(grid_param)
        except ValueError:
            raise WPSError('Error converting gaussian parameter to an int')

        grid = cdms2.createGaussianGrid(nlats)

        logger.info('Created target gaussian grid {}'.format(grid.shape))
    else:
        raise WPSError('Unknown grid type for regridding: {}', grid_type)

    return grid


def read_grid_from_file(gridder):
    try:
        with cdms2.open(gridder.grid.uri) as infile:
            data = infile(gridder.grid.var_name)
    except cdms2.CDMSError:
        raise WPSError('Failed to read the grid from {} in {}', gridder.grid.var_name, gridder.grid.uri)

    return data.getGrid()


def generate_grid(gridder, selector):
    if isinstance(gridder.grid, cwt.Variable):
        grid = read_grid_from_file(gridder)
    else:
        grid = generate_user_defined_grid(gridder)

    if isinstance(gridder.grid, cwt.Variable):
        grid_src = gridder.grid.uri
    else:
        grid_src = gridder.grid

    metrics.WPS_REGRID.labels(gridder.tool, gridder.method, grid_src).inc()

    return subset_grid(grid, selector)


def new_shape(shape, time_slice):
    diff = time_slice.stop - time_slice.start

    return (diff, ) + shape[1:]


def from_delayed(fm, uri, var, chunk_shape):
    size = var.shape[0]

    step = chunk_shape[0]

    logger.info('Building input from delayed functions size %r step %r', size, step)

    chunks = [slice(x, min(size, x+step), 1) for x in range(0, size, step)]

    logger.info('Generated %r chunks', len(chunks))

    delayed = [dask.delayed(retrieve_chunk)(uri, var.id, {'time': x}, fm.cert_data) for x in chunks]

    arrays = [da.from_delayed(x, new_shape(var.shape, y), var.dtype) for x, y in zip(delayed, chunks)]

    return da.concatenate(arrays, axis=0)


def slice_to_subaxis(axis_slice, axis):
    i = axis_slice.start or 0
    j = axis_slice.stop or len(axis)
    k = axis_slice.step or 1

    return i, j, k


class InputManager(object):
    def __init__(self, fm, uris, var_name):
        self.fm = fm
        self.uris = uris
        self.var_name = var_name
        self.attrs = {}
        self.vars = {}
        self.vars_axes = {}
        self.axes = {}
        self.target_var_axes = []
        self.units = None
        self.time_axis = []

    def __repr__(self):
        return ('InputManager(uris={!r}, var_name={!r}, attrs={!r}, '
                'vars={!r}, vars_axes={!r}, axes={!r})').format(self.uris, self.var_name,
                                                                self.attrs, self.vars,
                                                                self.vars_axes, self.axes)

    @classmethod
    def from_cwt_variable(cls, fm, variable):
        logger.info('Creating manager with single file')

        return cls(fm, [variable.uri, ], variable.var_name)

    @classmethod
    def from_cwt_variables(cls, fm, variables):
        var_names = [x.var_name for x in variables]

        if len(set(var_names)) > 1:
            raise WPSError('Mismatched variable names {!r}', var_names)

        uris = [x.uri for x in variables]

        logger.info('Creating manager with %r files', len(uris))

        return cls(fm, uris, var_names[0])

    @property
    def shape(self):
        return self.variable.shape

    @property
    def nbytes(self):
        return self.variable.nbytes

    @property
    def dtype(self):
        return self.variable.dtype

    @property
    def blocks(self):
        return self.variable.blocks

    @property
    def variable(self):
        return self.vars[self.var_name]

    @variable.setter
    def variable(self, value):
        self.vars[self.var_name] = value

    @property
    def variable_axes(self):
        return self.vars_axes[self.var_name]

    def remove_axis(self, axis_name):
        try:
            del self.axes[axis_name]
        except KeyError as e:
            raise WPSError('Did not find axis {!r}', e)
        else:
            for name in list(self.vars_axes.keys()):
                if name == self.var_name and axis_name in self.vars_axes[name]:
                    index = self.vars_axes[name].index(axis_name)

                    self.vars_axes[name].pop(index)
                else:
                    if axis_name in self.vars_axes[name]:
                        del self.vars_axes[name]

                    if name in self.vars:
                        del self.vars[name]

    def copy(self):
        new = InputManager(self.fm, self.uris, self.var_name)

        new.attrs = self.attrs.copy()

        for x, y in self.vars.items():
            new.vars[x] = y.copy()

        new.vars_axes = dict((x, y.copy()) for x, y in self.vars_axes.items())

        new.axes = dict((x, y.clone()) for x, y in self.axes.items())

        return new

    def load_axes(self, var, stored_time):
        axis_ids = []

        for axis in var.getAxisList():
            if self.var_name == var.id:
                self.target_var_axes.append(axis.id)

            if axis.isTime():
                if self.units is None:
                    self.units = axis.units

                    logger.info('Setting base units to %r', self.units)

                if not stored_time:
                    time_axis_clone = axis.clone()

                    time_axis_clone.toRelativeTime(self.units)

                    self.time_axis.append(time_axis_clone)

                    stored_time = True

                    self.attrs[axis.id] = axis.attributes.copy()

                    logger.info('Storing temporal axis %r shape %r', axis.id, axis.shape)
            elif axis.id not in self.axes:
                self.axes[axis.id] = axis.clone()

                self.attrs[axis.id] = axis.attributes.copy()

                logger.info('Storing spatial axis %r shape %r', axis.id, axis.shape)
            else:
                logger.info('Axis %r is not time and already discovered', axis.id)

            axis_ids.append(axis.id)

        logger.info('Setting variable %r axes to %r', var.id, axis_ids)

        self.vars_axes[var.id] = axis_ids

        return stored_time

    def load_variables(self, file, target_variable):
        stored_time = False

        for var in file.getVariables():
            # Skip all variables already discovered and do not contain temporal axis
            if var.id in self.vars and var.getTime() is None:
                logger.info('Skipping variable %r, already discovered')

                continue

            stored_time = self.load_axes(var, stored_time)

            self.attrs[var.id] = var.attributes.copy()

            if var.id == self.var_name:
                var_data = target_variable
            else:
                var_data = var()

            self.vars[var.id] = var_data

            logger.info('Storing variable %r shape %r', var.id, var.shape)

    def load_variables_and_axes(self, target_variable):
        for uri in self.uris:
            file = self.fm.open_file(uri)

            if 'global' not in self.attrs:
                logger.info('Storing global attributes')

                self.attrs['global'] = file.attributes.copy()

            self.load_variables(file, target_variable)

        if len(self.time_axis) > 1:
            self.axes['time'] = cdms2.MV2.axisConcatenate(self.time_axis, id='time', attributes=self.attrs['time'])

            logger.info('Updating temporal axis shape %r', self.axes['time'].shape)

            self.vars['time_bnds'] = cdms2.createVariable(
                self.axes['time'].getBounds(),
                axes=[self.axes[x] for x in self.vars_axes['time_bnds']],
                id='time_bnds',
                attributes=self.attrs['time_bnds'],
            )

            logger.info('Updating temporal bounds shape %r', self.vars['time_bnds'].shape)
        else:
            self.axes['time'] = self.time_axis[0]

    def subset_variables_and_axes(self, domain, target_variable):
        self.load_variables_and_axes(target_variable)

        map = map_domain(domain, self.axes)

        for name in list(self.axes.keys()):
            if name not in self.target_var_axes:
                continue

            axis_slice = map.get(name, slice(None, None, None))

            i, j, k = slice_to_subaxis(axis_slice, self.axes[name])

            shape = self.axes[name].shape

            self.axes[name] = self.axes[name].subAxis(i, j, k)

            logger.info('Subsetting axis %r shape %r -> %r', name, shape, self.axes[name].shape)

        for name in list(self.vars.keys()):
            selector = OrderedDict((x, map[x]) for x in self.vars_axes[name] if x in map)

            logger.info('Variable %r selector %r', name, selector)

            shape = self.vars[name].shape

            if name == self.var_name:
                dask_selector = tuple(selector.values())

                logger.info('Dask selector %r', dask_selector)

                self.vars[name] = self.vars[name][dask_selector]
            else:
                try:
                    self.vars[name] = self.vars[name](**selector)
                except TypeError:
                    # No Reason to subset
                    pass

            logger.info('Subsetting variable %r shape %r -> %r', name, shape, self.vars[name].shape)

        return self.vars[self.var_name]

    def to_xarray(self):
        axes = {}
        vars = {}

        logger.info('Building xarray with axes %r', self.axes.keys())
        logger.info('Building xarray with variables %r', self.vars.keys())

        for key, value in self.vars_axes.items():
            coords = OrderedDict()

            for axis_key in value:
                if axis_key in axes:
                    axis_value = axes[axis_key]
                else:
                    axis_value = axes[axis_key] = xr.DataArray(self.axes[axis_key], name=axis_key, dims=axis_key,
                                                               attrs=self.attrs[axis_key])

                    logger.info('Creating axis %r sizes %r', axis_value.name, axis_value.sizes)

                coords[axis_key] = axis_value

            vars[key] = xr.DataArray(self.vars[key], name=key, dims=coords.keys(), coords=coords, attrs=self.attrs[key])

            logger.info('Creating variable %r coords %r sizes %r', key, coords.keys(), vars[key].sizes)

        return xr.Dataset(vars, attrs=self.attrs['global'])

    def subset(self, domain):
        data = []

        logger.info('Subsetting the inputs')

        if len(self.uris) > 1:
            temporal_info = self.gather_temporal_info()

            self.sort_uris(temporal_info)

        for uri in self.uris:
            var = self.fm.get_variable(uri, self.var_name)

            chunks = (100, ) + var.shape[1:]

            if self.fm.requires_cert(uri):
                data.append(from_delayed(self.fm, uri, var, chunks))
            else:
                data.append(da.from_array(var, chunks=chunks))

            logger.info('Created input %r', data[-1])

        if len(data) > 1:
            data = da.concatenate(data, axis=0)

            logger.info('Concatenating data %r', data)
        else:
            data = data[0]

        subset_data = self.subset_variables_and_axes(domain, data)

        return data, subset_data

    def gather_temporal_info(self):
        ordering = []

        for uri in self.uris:
            var = self.fm.get_variable(uri, self.var_name)

            time = var.getTime()

            if time is None:
                raise WPSError('Unable to sort inputs {!r} has no time axis', uri)

            ordering.append((uri, time.units, time[0]))

        return ordering

    def sort_uris(self, ordering):
        logger.info('Sorting uris with %r', ordering)

        by_units = len(set([x[1] for x in ordering])) != 1

        logger.info('By Units: %r', by_units)

        by_first = len(set([x[2] for x in ordering])) != 1

        logger.info('By Firsts: %r', by_first)

        if by_units:
            ordering = sorted(ordering, key=lambda x: x[1])

            logger.info('Sorted uris by units')
        elif by_first:
            ordering = sorted(ordering, key=lambda x: x[2])

            logger.info('Sorted uris by first values')
        else:
            raise WPSError('Unable to determine ordering of files')

        self.uris = [x[0] for x in ordering]

        logger.info('Sorted uris %r', self.uris)
