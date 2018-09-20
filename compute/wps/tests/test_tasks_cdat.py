#! /usr/bin/env python

import os
import json
import datetime

import cwt
import mock
from django import test

from wps.tasks import cdat

# Need to use this to mock the return value of cdms2.open since 
class MockFile(object):
    def __init__(self, nbytes=0):
        self.called = False
        self.call_count = 0
        self.current = None
        self.returned = []
        self.mock_write = None
        self.nbytes = nbytes

    def write(self, *args, **kwargs):
        self.mock_write = (args, kwargs)

    def __call__(self, *args, **kwargs):
        self.called = True

        self.call_count += 1

        value = mock.MagicMock()

        type(value).nbytes = mock.PropertyMock(return_value=self.nbytes)

        self.returned.append((value, args, kwargs))

        return value

class CDATTaskTestCase(test.TestCase):

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.base.cdms2.MV2.concatenate')
    @mock.patch('os.stat')
    def test_concat_process_output(self, mock_stat, mock_mv, mock_open, mock_job):
        type(mock_stat.return_value).st_size = 1200000

        attrs = [
            {
                "file01.nc": {
                    'elapsed': 120,
                },
            },
            {
                "file02.nc": {
                    'elapsed': 60,
                },
            }
        ]

        input_paths = ['file01.nc', 'file02.nc']
    
        mock_file = MockFile()

        mock_open.return_value.__enter__.return_value = mock_file

        mock_op = mock.MagicMock()

        result = cdat.concat_process_output(attrs, input_paths, mock_op, 'tas', 'time', './output.nc', 0)

        mock_open.assert_any_call('file01.nc')
        mock_open.assert_any_call('file02.nc')

        self.assertEqual(mock_file.returned[0][1:], (('tas',), {}))
        self.assertEqual(mock_file.returned[1][1:], (('tas',), {}))

        mock_file.returned[0][0].getAxisIndex.assert_called_with('time')

        data = [x[0] for x in mock_file.returned]

        mock_mv.assert_called_with(data, axis=mock_file.returned[0][0].getAxisIndex.return_value)

        mock_open.assert_any_call('./output.nc', 'w')

        self.assertEqual(mock_file.mock_write, ((mock_mv.return_value,), {'id': 'tas'}))

        grouped = {}

        for x in attrs:
            grouped.update(x)

        self.assertEqual(result, grouped)

    @mock.patch('os.stat')
    def test_base_process_gridder(self, mock_stat):
        type(mock_stat.return_value).st_size = 1200000

        mapped = {
            'time': slice(0, 1),
            'lat': slice(20, 30),
            'lon': slice(10, 40),
        }

        attrs = {
            'file01.nc': {
                'path': './file01.nc',
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile(1e9)
        outfile = MockFile(1e9)

        now = datetime.datetime.now()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]
        mock_self.get_now.side_effect = [
            now,
            now+datetime.timedelta(minutes=30),
            now+datetime.timedelta(minutes=33),
            now+datetime.timedelta(minutes=40),
        ]

        op = cwt.Process('CDAT.subset')

        gridder = cwt.Gridder(grid='gaussian~32')

        op.parameters['gridder'] = gridder

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat'], './output.nc', 0)

        mock_self.generate_selector.assert_called_with(infile.returned[0][0])
        mock_self.generate_grid.assert_called_with(gridder)
        mock_self.subset_grid.assert_called_with(
            mock_self.generate_grid.return_value, 
            mock_self.generate_selector.return_value)
    
    def test_base_process_cached_gridder(self):
        mapped = {
            'time': slice(0, 1),
            'lat': slice(20, 30),
            'lon': slice(10, 40),
        }

        attrs = {
            'file01.nc': {
                'path': './file01.nc',
                'mapped': mapped,
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile(1e9)
        outfile = MockFile(1e9)

        now = datetime.datetime.now()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]
        mock_self.get_now.side_effect = [
            now,
            now+datetime.timedelta(minutes=30),
            now+datetime.timedelta(minutes=33),
            now+datetime.timedelta(minutes=40),
        ]

        op = cwt.Process('CDAT.subset')

        gridder = cwt.Gridder(grid='gaussian~32')

        op.parameters['gridder'] = gridder

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat'], './output.nc', 0)

        mock_self.generate_selector.assert_not_called()
        mock_self.generate_grid.assert_called_with(gridder)
        mock_self.subset_grid.assert_called_with(
            mock_self.generate_grid.return_value, mapped)

    def test_base_process_cached(self):
        mapped = {
            'time': slice(0, 1),
            'lat': slice(20, 30),
            'lon': slice(10, 40),
        }

        attrs = {
            'file01.nc': {
                'path': './file01.nc',
                'mapped': mapped,
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile(1e9)
        outfile = MockFile(1e9)

        now = datetime.datetime.now()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]
        mock_self.get_now.side_effect = [
            now,
            now+datetime.timedelta(minutes=30),
            now+datetime.timedelta(minutes=33),
            now+datetime.timedelta(minutes=40),
        ]

        op = cwt.Process('CDAT.subset')

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat'], './output.nc', 0)

        self.assertTrue(infile.called)
        self.assertEqual(infile.returned[0][1:], (('tas',), mapped))

    def test_base_process_average_multiple(self):
        attrs = {
            'file01.nc': {
                'path': './file01.nc',
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile()
        outfile = MockFile()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]

        op = cwt.Process('CDAT.average')

        op.add_parameters(weightoptions=['generate'])

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat', 'lon'], './output.nc', 0)

        axis_expected = ''.join([str(infile.returned[0][0].getAxisIndex.return_value) for _ in
                 range(2)])

        self.assertEqual(mock_self.PROCESS.call_count, 1)
        mock_self.PROCESS.assert_called_with(
            infile.returned[0][0], 
            axis=axis_expected,
            weights='generate')

    def test_base_process_average(self):
        attrs = {
            'file01.nc': {
                'path': './file01.nc',
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile()
        outfile = MockFile()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]

        op = cwt.Process('CDAT.average')

        op.add_parameters(weightoptions=['generate'])

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat'], './output.nc', 0)

        self.assertEqual(mock_self.PROCESS.call_count, 1)
        mock_self.PROCESS.assert_called_with(
            infile.returned[0][0], 
            axis=str(infile.returned[0][0].getAxisIndex.return_value),
            weights='generate')

    def test_base_process(self):
        attrs = {
            'file01.nc': {
                'path': './file01.nc',
            }
        }

        mock_self = mock.MagicMock()

        infile = MockFile(1e9)
        outfile = MockFile(1e9)

        now = datetime.datetime.now()

        mock_self.open.return_value.__enter__.side_effect = [infile, outfile]
        mock_self.get_now.side_effect = [
            now,
            now+datetime.timedelta(minutes=30),
            now+datetime.timedelta(minutes=33),
            now+datetime.timedelta(minutes=40),
        ]

        op = cwt.Process('CDAT.subset')

        cdat.base_process(mock_self, attrs, 'file01.nc', op, 'tas', 'days since 1990-01-01', ['lat'], './output.nc', 0)

        mock_self.load_job.assert_called_with(0)

        self.assertTrue(infile.called)
        self.assertEqual(infile.returned[0][1:], (('tas',), {}))

        infile.returned[0][0].regrid.assert_not_called()
        infile.returned[0][0].getAxisIndex.assert_called_with('lat')

        mock_self.generate_selector.assert_not_called()
        mock_self.generate_grid.assert_not_called()
        mock_self.subset_grid.assert_not_called()

        self.assertEqual(mock_self.PROCESS.call_count, 1)
        mock_self.PROCESS.assert_called_with(
            infile.returned[0][0], 
            axis=infile.returned[0][0].getAxisIndex.return_value)

        self.assertEqual(outfile.mock_write, ((mock_self.PROCESS.return_value,), {'id': 'tas'}))

    @mock.patch('wps.tasks.cdat.retrieve_data_cached')
    @mock.patch('os.stat')
    @mock.patch('os.makedirs')
    def test_base_retrieve_cached(self, mock_make, mock_stat, mock_retrieve):
        type(mock_stat.return_value).st_size = mock.PropertyMock(return_value=1200000)

        attrs = {
            'file01.nc': {
                'cached': True,
                'path': './file01.nc',
            }
        }

        now = datetime.datetime.now()

        mock_self = mock.MagicMock()

        mock_self.get_now.side_effect = [now,
                                         now+datetime.timedelta(seconds=30)]

        op = cwt.Process('CDAT.subset')

        cdat.base_retrieve(mock_self, attrs, ['file01.nc'], op, 'tas', 
                           'days since 1990-01-01', './output.nc', 0)

        mock_self.generate_grid.assert_not_called()
        mock_self.generate_selector.assert_not_called()
        mock_self.subset_grid.assert_not_called()

        self.assertEqual(mock_retrieve.call_count, 1)

    @mock.patch('wps.tasks.cdat.retrieve_data')
    @mock.patch('os.stat')
    @mock.patch('os.makedirs')
    def test_base_retrieve_regrid(self, mock_make, mock_stat, mock_retrieve):
        now = datetime.datetime.now()

        type(mock_stat.return_value).st_size = mock.PropertyMock(return_value=1200000)

        attrs = {
            'file01.nc': {
                'path': './file01.nc',
            }
        }

        mock_self = mock.MagicMock()

        mock_self.get_now.side_effect = [now,
                                         now+datetime.timedelta(seconds=30)]

        op = cwt.Process('CDAT.subset')

        gridder = cwt.Gridder(grid='gaussian~32')

        op.parameters['gridder'] = gridder

        cdat.base_retrieve(mock_self, attrs, ['file01.nc'], op, 'tas', 
                           'days since 1990-01-01', './output.nc', 0)

        mock_self.generate_grid.assert_called_with(gridder)

        mock_self.generate_selector.assert_called_with(
            mock_self.open.return_value.__enter__.return_value.__getitem__.return_value)

        mock_self.subset_grid.assert_called_with(
            mock_self.generate_grid.return_value, 
            mock_self.generate_selector.return_value)

    @mock.patch('wps.tasks.cdat.retrieve_data')
    @mock.patch('os.stat')
    @mock.patch('os.makedirs')
    def test_base_retrieve(self, mock_make, mock_stat, mock_retrieve):
        type(mock_stat.return_value).st_size = mock.PropertyMock(return_value=1200000)

        keys = ['file01.nc', 'file02.nc']
        
        attrs = [
            {
                keys[0]: {
                    'path': './file01.nc',
                },
            },
            {
                keys[1]: {
                    'path': './file02.nc',
                }
            }
        ]

        now = datetime.datetime.now()

        mock_self = mock.MagicMock()

        mock_self.get_now.side_effect = [now,
                                         now+datetime.timedelta(seconds=30)]

        op = cwt.Process('CDAT.subset')

        cdat.base_retrieve(mock_self, attrs, keys, op, 'tas', 'days since 1990-01-01', './output.nc', 0)

        self.assertEqual(mock_self.open.call_count, 3)        
        mock_self.open.assert_any_call('./output.nc', 'w')
        mock_self.open.assert_any_call('./file01.nc')
        mock_self.open.assert_any_call('./file02.nc')

        mock_self.generate_grid.assert_not_called()
        mock_self.generate_selector.assert_not_called()
        mock_self.subset_grid.assert_not_called()

        self.assertEqual(mock_retrieve.call_count, 2)

    @mock.patch('wps.tasks.cdat.retrieve_data')
    def test_retrieve_data_cached(self, mock_retrieve):
        now = datetime.datetime.now()

        type(mock_retrieve.return_value).nbytes = mock.PropertyMock(return_value=1e9)

        mock_self = mock.MagicMock()

        mock_self.get_now.side_effect = [
            now,
            now+datetime.timedelta(minutes=30),
            now+datetime.timedelta(minutes=32),
            now+datetime.timedelta(minutes=52),
        ]

        mock_infile = MockFile()

        mock_outfile = MockFile()

        mapped = {
            'time': slice(100, 200),
            'lat': slice(10, 24),
            'lon': slice(20, 30),
        }

        kwargs = {
            'mapped': mapped,
            'chunk_axis': 'time',
            'chunk_list': [slice(0, 1), slice(1, 2)]
        }

        grid = mock.MagicMock()

        gridder = cwt.Gridder(grid='gaussian~32', tool='ESMF', method='linear')

        cdat.retrieve_data_cached(mock_self, mock_infile, mock_outfile, 'tas', grid, gridder, 'days since 1990-01-01', **kwargs)

        mock_self.subset_grid.assert_called_with(grid, mapped)

    @mock.patch('wps.tasks.cdat.retrieve_data')
    def test_retrieve_data_cached(self, mock_retrieve):
        mock_self = mock.MagicMock()

        mock_infile = MockFile()

        mock_outfile = MockFile()

        mapped = {
            'time': slice(100, 200),
            'lat': slice(10, 24),
            'lon': slice(20, 30),
        }

        kwargs = {
            'mapped': mapped,
            'chunk_axis': 'time',
            'chunk_list': [slice(0, 1), slice(1, 2)]
        }

        cdat.retrieve_data_cached(mock_self, mock_infile, mock_outfile, 'tas', None, None, 'days since 1990-01-01', **kwargs)

        mock_self.subset_grid.assert_not_called()

        mock_retrieve.assert_called()
        self.assertEqual(mock_retrieve.call_count, 2)

        mapped.update({'time': slice(0, 1)})
        mock_retrieve.assert_any_call(mock_infile, mock_outfile, 'tas', None, None, 'days since 1990-01-01', mapped)

        mapped.update({'time': slice(1, 2)})
        mock_retrieve.assert_any_call(mock_infile, mock_outfile, 'tas', None, None, 'days since 1990-01-01', mapped)

    def test_retrieve_data_rebase_time(self):
        mock_infile = MockFile()

        mock_outfile = MockFile()

        cdat.retrieve_data(mock_infile, mock_outfile, 'tas', None, None, 'days since 1990-01-01')

        mock_infile.returned[0][0].getTime.assert_called()
        mock_infile.returned[0][0].getTime.return_value.toRelativeTime.assert_called_with('days since 1990-01-01')

    def test_retrieve_data_regrid(self):
        gridder = cwt.Gridder(grid='gaussian~32', tool='ESMF', method='linear')

        grid = mock.MagicMock()

        mock_infile = MockFile()

        mock_outfile = MockFile()

        cdat.retrieve_data(mock_infile, mock_outfile, 'tas', grid, gridder, None)

        mock_infile.returned[0][0].regrid.assert_called_with(grid, regridTool='ESMF', regridMethod='linear')

    def test_retrieve_data_selector(self):
        mock_infile = MockFile()

        mock_outfile = MockFile()

        mapped = {
            'time': slice(100, 200),
            'lat': slice(10, 24),
            'lon': slice(20, 30),
        }

        cdat.retrieve_data(mock_infile, mock_outfile, 'tas', None, None, None, mapped=mapped)

        self.assertEqual(mock_infile.returned[0][2], mapped)

    def test_retrieve_data(self):
        mock_infile = MockFile()

        mock_outfile = MockFile()

        cdat.retrieve_data(mock_infile, mock_outfile, 'tas', None, None, None)

        self.assertTrue(mock_infile.called)
        self.assertEqual(mock_infile.call_count, 1)
        self.assertEqual(mock_infile.returned[0][1], ('tas',))

        mock_infile.returned[0][0].regrid.assert_not_called()
        mock_infile.returned[0][0].getTime.assert_not_called()

        self.assertEqual(mock_outfile.mock_write[0], (mock_infile.returned[0][0],))
        self.assertEqual(mock_outfile.mock_write[1], {'id': 'tas'})
