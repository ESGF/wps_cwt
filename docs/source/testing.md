# Testing

* [Environment](#environment)
* [Running the tests](#running-the-tests)
  * [Django](#django)
  * [Angular](#angular)

##### Environment

Conda is the preferred testing environment. The following should all be
performed from the root

1. Install [Conda](https://conda.io/docs/user-guide/install/index.html)
2. Checkout source `git clone https://github.com/ESGF/esgf-compute-wps`
3. Move into the repo `cd esgf-compute-wps`
4. Create the testing environment `conda create env -n wps -f docker/common/environment.yaml`
5. Activate the environment `source activate wps`
6. Install the testing dependencies `pip install -r compute/wps/tests/requirements.txt`

##### Running the tests

###### Django

`python compute/manage.py test`

###### Angular

TODO
