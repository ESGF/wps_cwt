import mock
import requests
from django import test
from django.conf import settings

from wps import tasks
from wps import WPSError

class TaskMetricsTestCase(test.TestCase):

    @mock.patch.object(requests, 'get')
    @mock.patch('wps.models.Job.objects.filter')
    @mock.patch('wps.tasks.base.CWTBaseTask.load_job')
    def test_health(self, mock_job, mock_filter, mock_get):
        mock_get.return_value.json.return_value = {
            'status': 'ok', 
            'data': {
                'result': [
                    {
                        'metric': 'test',
                        'value': [],
                    },
                ]
            }
        }

        mock_filter.return_value.exclude.return_value.exclude.return_value.count.return_value = 2

        mock_filter.return_value.exclude.return_value.exclude.return_value.exclude.return_value.count.return_value = 2

        tasks.metrics(0, 1)

        mock_job.assert_called_with(1)

        self.assertEqual(mock_get.call_count, 15)

    @mock.patch.object(requests, 'get')
    def test_query_multiple_value_custom_type(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': [
                    {'metric': {'request': 'GetCapabilities'}, 'value': [32.2,
                                                                         '10.2'] },
                    {'metric': {'request': 'DescribeProcess'}, 'value': [32.2,
                                                                         '32.2'] },
                ]
            }
        }

        result = tasks.query_multiple_value('request', type=float, query='sum(wps_requests)')

        self.assertIsInstance(result, dict)
        self.assertEqual(result, {'GetCapabilities': 10.2, 'DescribeProcess':
                                  32.2})

    @mock.patch.object(requests, 'get')
    def test_query_multiple_value(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': [
                    {'metric': {'request': 'GetCapabilities'}, 'value': [32.2, 10] },
                    {'metric': {'request': 'DescribeProcess'}, 'value': [32.2, 32] },
                ]
            }
        }

        result = tasks.query_multiple_value('request', query='sum(wps_requests)')

        self.assertIsInstance(result, dict)
        self.assertEqual(result, {'GetCapabilities': 10, 'DescribeProcess':
                                  32})

    @mock.patch.object(requests, 'get')
    def test_query_single_value_custom_type(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': [
                    {'value': [32.2, 2.1]},
                ]
            }
        }

        result = tasks.query_single_value(type=float, query='sum(wps_requests)')

        self.assertIsInstance(result, float)
        self.assertEqual(result, 2.1)

    @mock.patch.object(requests, 'get')
    def test_query_single_value_empty_result(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': []
            }
        }

        result = tasks.query_single_value(query='sum(wps_requests)')

        self.assertIsInstance(result, int)
        self.assertEqual(result, 0)

    @mock.patch.object(requests, 'get')
    def test_query_single_value(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': [
                    {'value': [32.2, 2]},
                ]
            }
        }

        result = tasks.query_single_value(query='sum(wps_requests)')

        self.assertIsInstance(result, int)
        self.assertEqual(result, 2)

    @mock.patch.object(requests, 'get')
    def test_query_prometheus_malformed(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'data': {
                'result': {
                    'data': 3.2,
                }
            }
        }

        with self.assertRaises(WPSError):
            tasks.query_prometheus(query='sum(wps_requests)')

    @mock.patch.object(requests, 'get')
    def test_query_prometheus_not_ok(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=False)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': {
                    'data': 3.2,
                }
            }
        }

        with self.assertRaises(WPSError):
            tasks.query_prometheus(query='sum(wps_requests)')

    @mock.patch.object(requests, 'get')
    def test_query_prometheus(self, mock_get):
        type(mock_get.return_value).ok = mock.PropertyMock(return_value=True)

        mock_get.return_value.json.return_value = {
            'status': 'ok',
            'data': {
                'result': {
                    'data': 3.2,
                }
            }
        }

        result = tasks.query_prometheus(query='sum(wps_requests)')

        self.assertEqual(result, {'data': 3.2})

        mock_get.assert_called_with(settings.METRICS_HOST, params={'query':
                                                                   'sum(wps_requests)'})
