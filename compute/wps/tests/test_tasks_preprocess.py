#! /usr/bin/env python

import copy
import mock

import cwt
from django import test
from django.conf import settings

from wps import WPSError
from wps import helpers
from wps import tasks

class FakeCache(object):
    def __init__(self, items):
        self.items = items

    def count(self):
        return len(self.items)

    def __len__(self):
        return len(self.items)

    def __getitem__(self, index):
        return self.items[index]

class PreprocessTestCase(test.TestCase):

    def setUp(self):
        self.op = cwt.Process('CDAT.subset')
        self.op2 = cwt.Process('CDAT.average')

        lat = mock.MagicMock()
        lat.isTime.return_value = False
        lat.id = 'lat'
        lat.mapInterval.return_value = [100, 200]
        lat.__len__.return_value = 200
        lat.__getitem__.return_value = -90

        self.lat = lat

        lon = mock.MagicMock()
        lon.isTime.return_value = False
        lon.id = 'lon'
        lon.mapInterval.return_value = [50, 100]
        lon.__len__.return_value = 400
        lon.__getitem__.return_value = 0

        self.lon = lon

        time = mock.MagicMock()
        time.id = 'time'
        time.isTime.return_value = True
        time.clone.return_value.mapInterval.return_value = [200, 1000]
        time.__len__.return_value = 2000
        time.__getitem__.return_value = 30

        self.time = time

        time2 = mock.MagicMock()
        time2.id = 'time'
        time2.isTime.return_value = True
        time2.clone.return_value.mapInterval.return_value = [200, 1000]
        time2.__len__.return_value = 2000
        time2.__getitem__.return_value = 30

        self.time2 = time2

        time3 = mock.MagicMock()
        time3.id = 'time'
        time3.isTime.return_value = True
        time3.clone.return_value.mapIntervalue.return_value = [200, 1000]
        time3.__len__.return_value = 2000
        time3.__getitem__.return_value = 3600

        self.time3 = time3

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('requests.post')
    def test_wps_execute(self, mock_get, mock_job):
        attrs = [
            {
                'sort': 'units',
                'base_units': 'days since 1990-01-01',
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                    'cached': None,
                    'mapped': {
                        'lat': slice(100, 200, 1),
                        'lon': slice(50, 100, 1),
                        'time': slice(1000, 2000),
                    },
                    'chunks': {
                        'time': [slice(x, x+20) for x in xrange(1000, 2000, 20)],
                    }
                }
            },
            {
                'sort': 'units',
                'base_units': 'days since 1990-01-01',
                'file2.nc': {
                    'units': 'days since 1990-01-01',
                    'first': 30,
                    'cached': None,
                    'mapped': {
                        'lat': slice(100, 200, 1),
                        'lon': slice(50, 100, 1),
                        'time': slice(0, 1500),
                    },
                    'chunks': {
                        'time': [slice(x, x+20) for x in xrange(0, 1500, 20)],
                    }
                }
            }
        ]

        var1 = cwt.Variable('file1.nc', 'tas')
        var2 = cwt.Variable('file2.nc', 'tas')

        variable = { var1.name: var1, var2.name: var2 }

        d1 = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(1000, 3500))

        domain = { d1.name: d1 }

        op1 = cwt.Process('CDAT.subset')
        op1.inputs = [var1, var2]
        op1.domain = d1

        operation = { op1.name: op1 }

        tasks.wps_execute(attrs, variable, domain, operation, 0, 1)

        expected = {}

        for item in attrs:
            expected.update(item)

        expected['user_id'] = 0
        expected['job_id'] = 1
        expected['preprocess'] = True
        expected['root'] = op1.name
        expected['workflow'] = False
        expected['variable'] = variable
        expected['domain'] = domain
        expected['operation'] = operation

        expected = helpers.encoder(expected)

        mock_get.assert_called_with(settings.WPS_EXECUTE_URL, data=expected, verify=False)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_memory_check_average_multiple(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        # Should fail, max memory is 16000000 and axis dependency causes size
        # to be 100*50*1500*4=30000000, which is too large
        with self.settings(WORKER_MEMORY=1600000), self.assertRaises(WPSError):
            result = tasks.generate_chunks(attrs, self.op2, 'file1.nc', ['lat',
                                                                        'lon'], 0)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_memory_check_average(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        with self.settings(WORKER_MEMORY=1600000):
            result = tasks.generate_chunks(attrs, self.op2, 'file1.nc', ['lat'], 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'time': [slice(x, min(x+40, 1500), 1) for x in xrange(0, 1500, 40)],
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_memory_check_error(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        with self.settings(WORKER_MEMORY=2048), self.assertRaises(WPSError):
            tasks.generate_chunks(attrs, self.op, 'file1.nc', None, 0)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_cached(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': {
                    'path': './file1.nc',
                    'mapped': {
                        'lat': slice(20, 200, 1),
                        'lon': slice(10, 100, 1),
                        'time': slice(200, 1100, 1),
                    }
                },
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', None, 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'time': [slice(200, 1100, 1)],
        }

        self.assertEqual(result, expected)
        
    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_temporal(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', ['time'], 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'lat': [slice(100, 200, 1)],
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_spatial(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', ['lat'], 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'time': [slice(0, 1500, 1)],
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_temporal_spatial(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', ['time', 'lat'], 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'lon': [slice(50, 100, 1)],
        }

        self.assertEqual(result, expected)
        
    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_spatial_temporal(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', ['lat', 'time'], 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'lon': [slice(50, 100, 1)],
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks_multiple(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        with self.settings(WORKER_MEMORY=160000):
            result = tasks.generate_chunks(attrs, self.op, 'file1.nc', None, 0)

        print result

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {
            'time': [slice(x, x+4, 1) for x in xrange(0, 1500, 4)]
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_generate_chunks(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': None,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        result = tasks.generate_chunks(attrs, self.op, 'file1.nc', None, 0)

        expected = copy.deepcopy(attrs)

        expected['file1.nc']['chunks'] = {'time': [slice(0, 1500, 1)]}

        self.assertEqual(result, expected)

    @mock.patch('wps.models.Cache.objects.filter')
    def test_check_cache_entries_superset(self, mock_filter):
        mock_valid = mock.PropertyMock(return_value=True)

        mock_cache = mock.MagicMock()
        type(mock_cache).valid = mock_valid
        mock_cache.is_superset.return_value = True

        mock_filter.return_value = FakeCache([mock_cache])

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(200, 400))

        result = tasks.check_cache_entries('file1.nc', 'tas', domain, 0)

        self.assertEqual(result, mock_cache)

        mock_cache.is_superset.assert_called_with(domain)

        mock_cache.delete.assert_not_called()

        mock_valid.assert_called()

    @mock.patch('wps.models.Cache.objects.filter')
    def test_check_cache_entries_not_superset(self, mock_filter):
        mock_valid = mock.PropertyMock(return_value=True)

        mock_cache = mock.MagicMock()
        type(mock_cache).valid = mock_valid
        mock_cache.is_superset.return_value = False

        mock_filter.return_value = FakeCache([mock_cache])

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(200, 400))

        result = tasks.check_cache_entries('file1.nc', 'tas', domain, 0)

        self.assertEqual(result, None)

        mock_cache.delete.assert_not_called()

        mock_valid.assert_called()

    @mock.patch('wps.models.Cache.objects.filter')
    def test_check_cache_entries_not_valid(self, mock_filter):
        mock_valid = mock.PropertyMock(return_value=False)

        mock_cache = mock.MagicMock()
        type(mock_cache).valid = mock_valid

        mock_filter.return_value = FakeCache([mock_cache])

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(200, 400))

        result = tasks.check_cache_entries('file1.nc', 'tas', domain, 0)

        self.assertEqual(result, None)

        mock_cache.delete.assert_called()

        mock_valid.assert_called()

    def test_check_cache_entries(self):
        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(200, 400))

        result = tasks.check_cache_entries('file1.nc', 'tas', domain, 0)

        self.assertEqual(result, None)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.preprocess.check_cache_entries')
    def test_check_cache_matched(self, mock_entries, mock_job):
        data = {
            'lat': slice(20, 200),
            'lon': slice(40, 400),
            'time': slice(0, 2000),
        }

        self.maxDiff = None

        type(mock_entries.return_value).dimensions = mock.PropertyMock(return_value=helpers.encoder(data))
        type(mock_entries.return_value).local_path = mock.PropertyMock(return_value='./file1.nc')

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            },
        }

        result = tasks.check_cache(attrs, 'file1.nc', 'tas', 0)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'cached': {
                    'path': './file1.nc',
                    'mapped': {
                        'lat': slice(80, 180),
                        'lon': slice(10, 60),
                        'time': slice(0, 1500),
                    }
                },
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    def test_check_cache(self, mock_job):
        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            },
        }

        result = tasks.check_cache(attrs, 'file1.nc', 'tas', 0)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'cached': None,
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            }
        }

        self.assertEqual(result, expected)

    def test_map_domain_axis_map_error(self):
        self.lat.mapInterval.side_effect = TypeError()

        dimen = cwt.Dimension('lat', 100, 200, cwt.VALUES)

        result = tasks.map_domain_axis(self.lat, dimen, 'days since 1990-01-01')

        self.assertEqual(result, None)

    def test_map_domain_axis_indices(self):
        dimen = cwt.Dimension('lat', 100, 200, cwt.INDICES)

        result = tasks.map_domain_axis(self.lat, dimen, 'days since 1990-01-01')

        self.assertEqual(result, slice(100, 200, 1))

    def test_map_domain_axis_time(self):
        dimen = cwt.Dimension('time', 100, 200, cwt.VALUES)

        result = tasks.map_domain_axis(self.time, dimen, 'days since 1990-01-01')

        self.assertEqual(result, slice(200, 1000, 1))

    def test_map_domain_axis_unknown_crs(self):
        dimen = cwt.Dimension('lat', 100, 200, cwt.CRS('hello'))

        with self.assertRaises(WPSError):
            tasks.map_domain_axis(self.lat, dimen, 'days since 1990-01-01')

    def test_map_domain_axis(self):
        dimen = cwt.Dimension('lat', 100, 200, cwt.VALUES)

        result = tasks.map_domain_axis(self.lat, dimen, 'days since 1990-01-01')

        self.assertEqual(result, slice(100, 200, 1))

    def test_map_domain_axis_no_dimension(self):
        result = tasks.map_domain_axis(self.lat, None, 'days since 1990-01-01')

        self.assertEqual(result, slice(0, 200, 1))

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_map_domain_time_indices(self, mock_var, mock_open, mock_credentials, mock_job):
        self.maxDiff = None
        mock_var.return_value.getAxisListIndex.return_value = [0, 1, 2, 0, 1, 2]
        mock_var.return_value.getAxisIndex.side_effect = [0, 1, 2, 0, 1, 2]
        mock_var.return_value.getAxis.side_effect = [
            self.lat, self.lon, self.time,
            self.lat, self.lon, self.time2,
        ]

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
                'file2.nc': {
                    'units': 'days since 1990-01-01',
                    'first': 30,
                }
            }
        }

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=slice(1000, 3500))

        result = tasks.map_domain_time_indices(attrs, 'tas', domain, 0, 1)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(0, 1500, 1),
                }
            },
            'file2.nc': {
                'units': 'days since 1990-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(1000, 2000, 1),
                }
            }
        }

        self.assertEqual(result, expected)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_map_domain_none(self, mock_var, mock_open, mock_credentials, mock_job):
        mock_var.return_value.getAxisList.return_value = [self.lat, self.lon, self.time]

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
            }
        }

        result = tasks.map_domain(attrs, 'file1.nc', 'tas', None, 0, 1)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(0, 200, 1),
                    'lon': slice(0, 400, 1),
                    'time': slice(0, 2000, 1),
                }
            }
        }

        self.assertEqual(result, expected)

        mock_open.assert_called_with('file1.nc')
        
        # Can't verify first called arg because of how mock handles __enter__
        self.assertEqual(mock_var.call_args[0][1], 'tas')

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_map_domain_missing_dimension(self, mock_var, mock_open, mock_credentials, mock_job):
        mock_var.return_value.getAxisIndex.side_effect = [0, -1, 2]
        mock_var.return_value.getAxis.side_effect = [self.lat, self.lon, self.time]

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
            }
        }

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(400, 800), name='d0')

        with self.assertRaises(WPSError):
            tasks.map_domain(attrs, 'file1.nc', 'tas', domain, 0, 1)

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_map_domain_partial(self, mock_var, mock_open, mock_credentials, mock_job):
        mock_var.return_value.getAxisListIndex.return_value = [0, 1, 2]
        mock_var.return_value.getAxisIndex.side_effect = [1, 2]
        mock_var.return_value.getAxis.side_effect = [self.lon, self.time, self.lat]

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
            }
        }

        domain = cwt.Domain(lon=(90, 270), time=(400, 800), name='d0')

        result = tasks.map_domain(attrs, 'file1.nc', 'tas', domain, 0, 1)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(0, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(200, 1000, 1),
                }
            }
        }

        self.assertEqual(result, expected)

        mock_open.assert_called_with('file1.nc')
        
        # Can't verify first called arg because of how mock handles __enter__
        self.assertEqual(mock_var.call_args[0][1], 'tas')

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_map_domain(self, mock_var, mock_open, mock_credentials, mock_job):
        mock_var.return_value.getAxisListIndex.return_value = [0, 1, 2]
        mock_var.return_value.getAxisIndex.side_effect = [0, 1, 2]
        mock_var.return_value.getAxis.side_effect = [self.lat, self.lon, self.time]

        attrs = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
            }
        }

        domain = cwt.Domain(lat=(-45, 45), lon=(90, 270), time=(400, 800), name='d0')

        result = tasks.map_domain(attrs, 'file1.nc', 'tas', domain, 0, 1)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'file1.nc': {
                'units': 'days since 2001-01-01',
                'first': 30,
                'mapped': {
                    'lat': slice(100, 200, 1),
                    'lon': slice(50, 100, 1),
                    'time': slice(200, 1000, 1),
                }
            }
        }

        self.assertEqual(result, expected)

        mock_open.assert_called_with('file1.nc')
        
        # Can't verify first called arg because of how mock handles __enter__
        self.assertEqual(mock_var.call_args[0][1], 'tas')

    @mock.patch('wps.tasks.CWTBaseTask.load_job')
    @mock.patch('wps.tasks.CWTBaseTask.load_credentials')
    @mock.patch('wps.tasks.CWTBaseTask.open')
    @mock.patch('wps.tasks.CWTBaseTask.get_variable')
    def test_determine_base_units(self, mock_var, mock_open, mock_credentials, mock_job):
        type(mock_var.return_value.getTime.return_value).units = mock.PropertyMock(side_effect=[
            'days since 2001-01-01',
            'days since 1990-01-01',
        ])

        mock_var.return_value.getTime.return_value.__getitem__.side_effect = [
            30, 30
        ]

        urls = ['file1.nc', 'file2.nc']

        result = tasks.determine_base_units(urls, 'tas', 0, 1)

        mock_job.assert_called_with(1)

        mock_credentials.assert_called_with(0)

        expected = {
            'sort': 'units',
            'base_units': 'days since 1990-01-01',
            'files': {
                'file1.nc': {
                    'units': 'days since 2001-01-01',
                    'first': 30,
                }, 
                'file2.nc': {
                    'units': 'days since 1990-01-01',
                    'first': 30,
                }
            }
        }

        self.assertEqual(result, expected)
