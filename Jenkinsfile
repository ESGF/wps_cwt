pipeline {
  agent none
  stages {
    stage('Build') {
      parallel {
        stage('provisioner') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              expression {
                return params.FORCE_PROVISIONER
              }

              changeset '**/compute_provisioner/**'
            }

          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make provisioner REGISTRY=${OUTPUT_REGISTRY}'
            }

          }
        }

        stage('tasks') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              expression {
                return params.FORCE_TASKS
              }

              changeset '**/compute_tasks/**'
            }

          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make tasks REGISTRY=${OUTPUT_REGISTRY}'
            }

          }
        }

        stage('wps') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              expression {
                return params.FORCE_WPS
              }

              changeset '**/compute_wps/**'
            }

          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make wps REGISTRY=${OUTPUT_REGISTRY}'
            }

          }
        }

        stage('thredds') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              expression {
                return params.FORCE_THREDDS
              }

              changeset '**/docker/thredds/**'
            }

          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make thredds REGISTRY=${OUTPUT_REGISTRY}'
            }

          }
        }

      }
    }

    stage('Testing') {
      parallel {
        stage('tasks') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              changeset '**/compute_tasks/**'
            }

          }
          environment {
            MPC = credentials('myproxyclient')
          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make tasks TARGET=testresult'
              sh 'chown -R 10000:10000 output'
            }

            junit(testResults: 'output/unittesting.xml', healthScaleFactor: 1)
            cobertura(coberturaReportFile: 'output/coverage.xml', autoUpdateStability: true, autoUpdateHealth: true, failNoReports: true, failUnhealthy: true, failUnstable: true, maxNumberOfBuilds: 2)
          }
        }

        stage('wps') {
          agent {
            node {
              label 'jenkins-buildkit'
            }

          }
          when {
            anyOf {
              changeset '**/compute_wps/**'
            }

          }
          steps {
            container(name: 'buildkit', shell: '/bin/sh') {
              sh 'make wps TARGET=testresult'
              sh 'chown -R 10000:10000 output'
            }

            junit(testResults: 'output/unittest.xml', healthScaleFactor: 1)
            cobertura(coberturaReportFile: 'output/coverage.xml', autoUpdateHealth: true, maxNumberOfBuilds: 2, failNoReports: true, failUnhealthy: true, failUnstable: true, autoUpdateStability: true)
          }
        }

      }
    }

    stage('Deploy') {
      agent {
        node {
          label 'jenkins-helm'
        }

      }
      when {
        anyOf {
          branch 'master'
          branch 'devel'
        }

      }
      environment {
        GH = credentials('ae3dd8dc-817a-409b-90b9-6459fb524afc')
      }
      steps {
        container(name: 'helm', shell: '/bin/bash') {
          sh '''#! /bin/bash

GIT_DIFF="$(git diff --name-only ${GIT_COMMIT} ${GIT_PREVIOUS_COMMIT})"

echo -e "GIT_DIFF\\n${GIT_DIFF}"

git clone -b devel https://github.com/esgf-compute/charts

cd charts/

conda install -c conda-forge ruamel.yaml

if [[ ! -z "$(echo ${GIT_DIFF} | grep /compute_provisioner/)" ]] || [[ "${FORCE_PROVISIONER}" == "true" ]]
then
  python scripts/update_config.py compute/values.yaml provisioner ${GIT_COMMIT:0:8}
fi

if [[ ! -z "$(echo ${GIT_DIFF} | grep /compute_wps/)" ]] || [[ "${FORCE_WPS}" == "true" ]]
then
  python scripts/update_config.py compute/values.yaml wps ${GIT_COMMIT:0:8}
fi

if [[ ! -z "$(echo ${GIT_DIFF} | grep /compute_tasks/)" ]] || [[ "${FORCE_TASKS}" == "true" ]]
then
  python scripts/update_config.py compute/values.yaml celery ${GIT_COMMIT:0:8}
fi

if [[ ! -z "$(echo ${GIT_DIFF} | grep /docker/thredds/)" ]] || [[ "${FORCE_THREDDS}" == "true" ]]
then
  python scripts/update_config.py compute/values.yaml thredds ${GIT_COMMIT:0:8}
fi

git config user.email ${GIT_EMAIL}

git config user.name ${GIT_NAME}

git add compute/values.yaml

git commit -m "Updates imageTag to ${GIT_COMMIT:0:8}"

git push https://${GH_USR}:${GH_PSW}@github.com/esgf-compute/charts'''
        }

      }
    }

  }
  parameters {
    booleanParam(name: 'FORCE_PROVISIONER', defaultValue: false, description: 'Force provisioner build')
    booleanParam(name: 'FORCE_TASKS', defaultValue: false, description: 'Force tasks(celery) build')
    booleanParam(name: 'FORCE_WPS', defaultValue: false, description: 'Force WPS build')
    booleanParam(name: 'FORCE_THREDDS', defaultValue: false, description: 'Force THREDDS build')
  }
}