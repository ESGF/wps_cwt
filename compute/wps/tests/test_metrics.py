from django import test

import mock
from wps import metrics

import cwt

class MetricsTestCase(test.TestCase):
    @mock.patch('wsgiref.simple_server.make_server')
    @mock.patch('prometheus_client.make_wsgi_app')
    def test_serve_metrics(self, mock_app, mock_make):
        metrics.serve_metrics()

        mock_make.assert_called_with('0.0.0.0', 8080, mock_app.return_value)

    @mock.patch('wps.metrics.inspect')
    def test_jobs_running(self, mock_inspect):
        mock_inspect.return_value.active.side_effect = AttributeError()

        result = metrics.jobs_running()

        self.assertEqual(result, 0)

    @mock.patch('wps.metrics.inspect')
    def test_jobs_running(self, mock_inspect):
        mock_inspect.return_value.active.return_value.values.return_value = [[{}]]

        result = metrics.jobs_running()

        self.assertEqual(result, 1)

    @mock.patch('wps.metrics.inspect')
    def test_jobs_queued_empty(self, mock_inspect):
        mock_inspect.return_value.scheduled.side_effect = AttributeError()

        result = metrics.jobs_queued()

        self.assertEqual(result, 0)

    @mock.patch('wps.metrics.inspect')
    def test_jobs_queued(self, mock_inspect):
        mock_inspect.return_value.scheduled.return_value.values.return_value = [[{}]]

        mock_inspect.return_value.reserved.return_value.values.return_value = [[{}]]

        result = metrics.jobs_queued()

        self.assertEqual(result, 2)

    @mock.patch('wps.metrics.FILE_ACCESSED')
    def test_track_file(self, mock_metric):
        var = cwt.Variable('http://somesite.com/doesntexist.nc', 'tas')

        metrics.track_file(var)

        mock_metric.labels.assert_called_with('somesite.com',
                                              '/doesntexist.nc',
                                              'http://somesite.com/doesntexist.nc',
                                              'tas')

        mock_metric.labels.return_value.inc.assert_called_once()
