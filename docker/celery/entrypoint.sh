#! /bin/bash

source activate wps

pushd $CWT_BASE

exec celery worker -A compute ${@} &

celery_ret=${?}
celery_pid=${!}

if [[ ${celery_ret} -ne 0 ]]; then
  echo "Celery failed to start"

  return ${celery_ret}
fi

if [[ -n "${CWT_METRICS}" ]]; then
  [[ ! -e "${CWT_METRICS}" ]] && mkdir "${CWT_METRICS}"

  exec python wps/metrics.py &

  metrics_ret=${?}
  metrics_pid=${!}

  if [[ ${metrics_ret} -ne 0 ]]; then
    echo "Metrics failed to start"

    return ${metrics_ret}
  fi
fi

trap "{ kill -9 ${celery_pid}; kill -9 ${metrics_pid}; rm -f ${CWT_METRICS}/*gauge*;" SIGINT SIGTERM

while sleep 60; do
  if [[ $(kill -0 ${celery_pid}) -ne 0 ]]; then
    echo "Celery process died"

    exit 1
  fi

  if [[ -n "${CWT_METRICS}" ]] && [[ $(kill -0 ${metrics_pid}) -ne 0 ]]; then
    echo "Metrics process died"

    exit 1
  fi
done
