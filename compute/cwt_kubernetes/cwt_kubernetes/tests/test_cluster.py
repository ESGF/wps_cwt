import mock
import unittest

from kubernetes import client

from cwt_kubernetes.cluster import Cluster


class TestCluster(unittest.TestCase):

    @mock.patch('cwt_kubernetes.cluster.config')
    def test_scale_down_exception(self, config):
        mock_pod = mock.MagicMock()

        cluster = Cluster('default', mock_pod)

        namespaced_pods = mock.MagicMock()

        mock_pod1 = mock.MagicMock()
        mock_pod1.status.pod_ip = '192.168.0.2'
        mock_pod1.metadata.name = 'dask-worker-test2'

        mock_pod2 = mock.MagicMock()
        mock_pod2.status.pod_ip = '192.168.0.3'
        mock_pod2.metadata.name = 'dask-worker-test3'

        mock_pod3 = mock.MagicMock()
        mock_pod3.status.pod_ip = '192.168.0.4'
        mock_pod3.metadata.name = 'dask-worker-test4'

        type(namespaced_pods).items = mock.PropertyMock(return_value=[
            mock_pod1,
            mock_pod2,
            mock_pod3,
        ])

        with mock.patch.object(cluster, 'core') as mock_core:
            mock_core.list_namespaced_pod.return_value = namespaced_pods

            mock_core.delete_namespaced_pod.side_effect = [
                client.rest.ApiException,
                None
            ]

            cluster.scale_down(['192.168.0.4', '192.168.0.2'])

            mock_core.list_namespaced_pod.assert_called_with('default')

            mock_core.delete_namespaced_pod.assert_any_call('dask-worker-test2', 'default')

            mock_core.delete_namespaced_pod.assert_any_call('dask-worker-test4', 'default')

    @mock.patch('cwt_kubernetes.cluster.config')
    def test_scale_down_empty(self, config):
        mock_pod = mock.MagicMock()

        cluster = Cluster('default', mock_pod)

        with mock.patch.object(cluster, 'core') as mock_core:
            cluster.scale_down([])

            mock_core.delete_namespaced_pod.assert_not_called()

    @mock.patch('cwt_kubernetes.cluster.config')
    def test_scale_down(self, config):
        mock_pod = mock.MagicMock()

        cluster = Cluster('default', mock_pod)

        namespaced_pods = mock.MagicMock()

        mock_pod1 = mock.MagicMock()
        mock_pod1.status.pod_ip = '192.168.0.2'
        mock_pod1.metadata.name = 'dask-worker-test2'

        mock_pod2 = mock.MagicMock()
        mock_pod2.stauts.pod_ip = '192.168.0.3'
        mock_pod2.metadata.name = 'dask-worker-test3'

        type(namespaced_pods).items = mock.PropertyMock(return_value=[
            mock_pod1,
            mock_pod2,
        ])

        with mock.patch.object(cluster, 'core') as mock_core:
            mock_core.list_namespaced_pod.return_value = namespaced_pods

            cluster.scale_down(['192.168.0.1', '192.168.0.2'])

            mock_core.list_namespaced_pod.assert_called_with('default')

            mock_core.delete_namespaced_pod.assert_called_with('dask-worker-test2', 'default')

    @mock.patch('cwt_kubernetes.cluster.config')
    @mock.patch('uuid.uuid4')
    def test_scale_up_exception(self, mock_uuid4, config):
        mock_pod = mock.MagicMock()

        mock_uuid4.side_effect = ['test1', 'test2']

        cluster = Cluster('default', mock_pod)

        with mock.patch.object(cluster, 'core') as mock_core:
            mock_core.create_namespaced_pod.side_effect = client.rest.ApiException()

            with self.assertRaises(client.rest.ApiException):
                cluster.scale_up(2)

    @mock.patch('cwt_kubernetes.cluster.config')
    @mock.patch('uuid.uuid4')
    def test_scale_up(self, mock_uuid4, config):
        mock_pod = mock.MagicMock()

        mock_uuid4.side_effect = ['test1', 'test2']

        cluster = Cluster('default', mock_pod)

        with mock.patch.object(cluster, 'core') as mock_core:
            cluster.scale_up(2)

            mock_core.create_namespaced_pod.assert_any_call('default', mock_pod)

    @mock.patch('cwt_kubernetes.cluster.config')
    def test_from_yaml(self, config):
        cluster = Cluster.from_yaml('default', 'tests/worker-spec.yml')

        self.assertIsNotNone(cluster.pod)