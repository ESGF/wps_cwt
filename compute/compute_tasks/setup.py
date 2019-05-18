import os
import setuptools

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setuptools.setup(
    name='compute-tasks',
    version='devel',
    author='Jason Boutte',
    author_email='boutte3@llnl.gov',
    description='Celery and Dask compute tasks',
    url='https://github.com/ESGF/esgf-compute-wps',
    packages=['compute_tasks'],
    entry_points={
        'console_scripts': [
            'compute-tasks-metrics=compute_tasks.metrics_:main',
            'compute-tasks-backend=compute_tasks.backend:main',
        ],
    }
)
