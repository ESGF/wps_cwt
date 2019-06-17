import os

FRONTEND_PORT = os.environ.get('FRONTEND_PORT', 7777)
BACKEND_PORT = os.environ.get('BACKEND_PORT', 7778)
REDIS_HOST = os.environ.get('REDIS_HOST', '')

ENCODING = 'ISO-8859-1'

HEARTBEAT_LIVENESS = 3
HEARTBEAT_INTERVAL = 1.0

INTERVAL_INIT = 1
INTERVAL_MAX = 32

READY = b'READY'
HEARTBEAT = b'HEARTBEAT'
REQUEST = b'REQUEST'
RESOURCE = b'RESOURCE'
ACK = b'ACK'
ERR = b'ERR'

ACK_TIMEOUT = 8
