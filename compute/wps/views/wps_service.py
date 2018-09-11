import json
import urllib
import re
import StringIO

import celery
import cwt
import django
from django import http
from django.conf import settings
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.decorators.http import require_http_methods
from lxml import etree

import wps
from . import common
from wps import backends
from wps import helpers
from wps import metrics
from wps import models
from wps import tasks
from wps import WPSError

logger = common.logger

class WPSExceptionError(WPSError):
    def __init__(self, message, code):
        self.report = wps.exception_report(message, code)

        super(WPSExceptionError, self).__init__('WPS Exception')

class WPSScriptGenerator(object):
    def __init__(self, variable, domain, operation, user):
        self.operation = operation

        self.domain = domain

        self.variable = variable

        self.user = user

        self.root_op = None

    def write_header(self, buf):
        buf.write('import cwt\nimport time\n\n')

        buf.write('api_key=\'{}\'\n\n'.format(self.user.auth.api_key))

        buf.write('wps = cwt.WPSClient(\'{}\', api_key=api_key)\n\n'.format(settings.WPS_ENDPOINT))

    def write_variables(self, buf):
        for v in self.variable.values():
            buf.write('var-{} = cwt.Variable(\'{}\', \'{}\')\n'.format(v.name, v.uri, v.var_name))

        buf.write('\n')

    def write_domain(self, buf):
        for d in self.domain.values():
            buf.write('dom-{} = cwt.Domain([\n'.format(d.name))

            for dim in d.dimensions:
                buf.write('\tcwt.Dimension(\'{}\', \'{}\', \'{}\', crs=cwt.VALUES, step=\'{}\'),\n'.format(dim.name, dim.start, dim.end, dim.step))

            buf.write('])\n\n')

    def write_processes(self, buf):
        for o in self.operation.values():
            buf.write('op-{} = wps.process(\'{}\')\n\n'.format(o.name, o.identifier))

        op = self.operation.values()[0]

        self.root_op = 'op-{}'.format(op.name)

        buf.write('wps.execute(op-{}, inputs=['.format(op.name))

        l = len(op.inputs) - 1

        for i, v in enumerate(op.inputs):
            buf.write('var-{}'.format(v.name))

            if i < l:
                buf.write(', ')

        if op.domain is not None:
            buf.write('], domain=dom-{}'.format(op.domain.name))

        if 'domain' in op.parameters:
            del op.parameters['domain']

        if 'gridder' in op.parameters:
            g = op.parameters['gridder']

            buf.write(', gridder=cwt.Gridder(\'{}\', \'{}\', \'{}\')'.format(g.tool, g.method, g.grid))

            del op.parameters['gridder']

        for p in op.parameters.values():
            buf.write(', {}={}'.format(p.name, p.values))

        buf.write(')\n\n')

    def write_status(self, buf):
        buf.write('while {}.processing:\n'.format(self.root_op))

        buf.write('\tprint {}.status\n\n'.format(self.root_op))

        buf.write('\ttime.sleep(1)\n\n'.format(self.root_op))

    def write_output(self, buf):
        buf.write('print {}.status\n\n'.format(self.root_op))

        buf.write('try:\n\timport vcs\n\timport cdms2\nexcept:\n\tpass\nelse:\n')

        buf.write('\tf = cdms2.open({}.output.uri)\n\n'.format(self.root_op))

        buf.write('\tv = vcs.init()\n\n')

        buf.write('\tv.plot(f[{}.output.var_name])'.format(self.root_op))

    def generate(self):
        buf = StringIO.StringIO()

        self.write_header(buf)

        self.write_variables(buf)

        self.write_domain(buf)

        self.write_processes(buf)

        self.write_status(buf)

        self.write_output(buf)

        data = buf.getvalue()

        buf.close()

        return data

def load_data_inputs(data_inputs, resolve_inputs=False):
    o, d, v = cwt.WPSClient.parse_data_inputs(data_inputs)

    v = dict((x.name, x) for x in v)

    d = dict((x.name, x) for x in d)

    o = dict((x.name, x) for x in o)

    logger.info('Loaded variables %r', v)

    logger.info('Loaded domains %r', d)

    logger.info('Loaded operations %r', o)

    if resolve_inputs:
        collected_inputs = list(y for x in o.values() for y in x.inputs)

        try:
            root_op = [o[x] for x in o.keys() if x not in collected_inputs][0]
        except IndexError as e:
            raise WPSError('Error resolving root operation')

        root_op.resolve_inputs(v, o)

        try:
            for x in o.values():
                if x.domain is not None:
                    x.domain = d[x.domain]
        except KeyError as e:
            raise WPSerror('Error resolving domain')

    return o, d, v

def get_parameter(params, name, required=True):
    """ Gets a parameter from a django QueryDict """

    # Case insesitive
    temp = dict((x.lower(), y) for x, y in params.iteritems())

    if name.lower() not in temp:
        logger.info('Missing required parameter %s', name)

        raise WPSExceptionError(name.lower(), cwt.ows.MissingParameterValue)

    param = temp.get(name.lower())

    if required and param is None:
        raise WPSError('Missing required parameter')

    return param

def handle_get_capabilities(host=None):
    if host is None:
        host = 'default'

    server = models.Server.objects.get(host=host)

    return server.capabilities

def handle_describe_process(identifiers):
    processes = models.Process.objects.filter(identifier__in=identifiers)

    return wps.process_descriptions_from_processes(processes)

def handle_execute(api_key, identifier, data_inputs):
    try:
        user = models.User.objects.filter(auth__api_key=api_key)[0]
    except IndexError:
        raise WPSError('Missing API key for WPS execute request')

    try:
        process = models.Process.objects.get(identifier=identifier)
    except models.Process.DoesNotExist:
        raise WPSError('Process "{identifier}" does not exist', identifier=identifier)

    # load the process drescription to get the data input descriptions
    process_descriptions = cwt.wps.CreateFromDocument(process.description)

    description = process_descriptions.ProcessDescription[0]

    kwargs = {}

    # load up the required data inputs for the process
    if description.DataInputs is not None:
        for x in description.DataInputs.Input:
            input_id = x.Identifier.value().lower()

            try:
                kwargs[input_id] = data_inputs[input_id]
            except KeyError:
                raise WPSError('Missing required input "{input_id}" for process {name}', input_id=input_id, name=identifier)

    server = models.Server.objects.get(host='default')

    job = models.Job.objects.create(server=server, process=process, user=user, extra=json.dumps(data_inputs))

    kwargs.update({
        'identifier': identifier,
        'user': user,
        'job': job,
        'process': process,
    })

    try:
        backend = backends.Backend.get_backend(process.backend)

        if backend is None:
            raise WPSError('Unknown backend "{name}"', name=process.backend)

        backend.execute(**kwargs)
    except Exception as e:
        job.failed(str(e))

        raise

    # at this point we've accepted the job
    job.accepted()

    return job.report

def handle_get(params, query_string):
    """ Handle an HTTP GET request. """
    request = get_parameter(params, 'request', True).lower()

    service = get_parameter(params, 'service', True)

    data_inputs = {}

    logger.info('GET request %r, service %r', request, service)

    if request == 'getcapabilities':
        with metrics.WPS_REQUESTS.labels('GetCapabilities', 'GET').time():
            response = handle_get_capabilities() 
    elif request == 'describeprocess':
        identifier = get_parameter(params, 'identifier', True).split(',')

        with metrics.WPS_REQUESTS.labels('DescribeProcess', 'GET').time():
            response = handle_describe_process(identifier)
    elif request == 'execute':
        api_key = get_parameter(params, 'api_key', True)

        identifier = get_parameter(params, 'identifier', True)

        query_string = urllib.unquote(query_string)

        if 'datainputs' in query_string:
            match = re.match('.*datainputs=\[(.*)\].*', query_string)

            if match is not None:
                # angular2 encodes ; breaking django query_string parsing so the 
                # webapp replaces ; with | and the change is reverted before parsing
                # the datainputs
                data_inputs = re.sub('\|(operation|domain|variable)=', ';\\1=', match.group(1))

                data_inputs = dict(x.split('=') for x in data_inputs.split(';'))

        with metrics.WPS_REQUESTS.labels('Execute', 'GET').time():
            response = handle_execute(api_key, identifier, data_inputs)
    else:
        raise WPSError('Operation "{name}" is not supported', name=request)

    return response

def handle_post(data, params):
    """ Handle an HTTP POST request. 

    NOTE: we only support execute requests as POST for the moment
    """
    try:
        request = cwt.wps.CreateFromDocument(data)
    except Exception as e:
        raise WPSError('Malformed WPS execute request')

    if isinstance(request, cwt.wps.CTD_ANON_11):
        raise WPSError('GetCapabilities POST not supported')
    elif isinstance(request, cwt.wps.CTD_ANON_12):
        raise WPSError('DescribeProcess POST not supported')

    data_inputs = {}

    for x in request.DataInputs.Input:
        data_inputs[x.Identifier.value()] = x.Data.LiteralData.value()

    api_key = params.get('api_key')

    logger.info('Handling POST request for API key %s', api_key)

    with metrics.WPS_REQUESTS.labels('Execute', 'POST').time():
        response = handle_execute(api_key, request.Identifier.value(), data_inputs)

    return response

@metrics.WPS_ERRORS.count_exceptions()
def handle_request(request):
    """ Convert HTTP request to intermediate format. """
    if request.method == 'GET':
        return handle_get(request.GET, request.META['QUERY_STRING'])
    elif request.method == 'POST':
        return handle_post(request.body, request.GET)

@require_http_methods(['POST'])
def generate(request):
    try:
        common.authentication_required(request)

        data_inputs = request.POST['datainputs']

        data_inputs = re.sub('\|(domain|operation|variable)=', ';\\1=', data_inputs, 3)

        o, d, v = load_data_inputs(data_inputs, resolve_inputs=True)

        script = WPSScriptGenerator(v, d, o, request.user)

        kernel = o.values()[0].identifier.split('.')[1]

        data = {
            'filename': '{}.py'.format(kernel),
            'text': script.generate()
        }
    except WPSError as e:
        return common.failed(str(e))
    else:
        return common.success(data)

@require_http_methods(['GET', 'POST'])
@ensure_csrf_cookie
def wps_entrypoint(request):
    response = None

    try:
        response = handle_request(request)
    except WPSExceptionError as e:
        logger.exception('WSPExceptionError')

        response = e.report
    except WPSError as e:
        logger.exception('WPSError')

        response = wps.exception_report(str(e), cwt.ows.NoApplicableCode)
    except Exception as e:
        logger.exception('Some generic exception')

        error = 'Please copy the error and report on Github: {}'.format(str(e))

        response = wps.exception_report(error, cwt.ows.NoApplicableCode)

    return http.HttpResponse(response, content_type='text/xml')

@require_http_methods(['GET'])
def regen_capabilities(request):
    try:
        common.authentication_required(request)

        common.authorization_required(request)

        server = models.Server.objects.get(host='default')

        processes = []

        for process in server.processes.all():
            proc = cwt.wps.process(process.identifier, process.identifier, '1.0.0')

            processes.append(proc)

        process_offerings = cwt.wps.process_offerings(processes)

        server.capabilities = wps.generate_capabilities(process_offerings)

        server.save()
    except WPSError as e:
        return common.failed(str(e))
    else:
        return common.success('Regenerated capabilities')

@require_http_methods(['GET'])
def status(request, job_id):
    try:
        job = models.Job.objects.get(pk=job_id)
    except models.Job.DoesNotExist:
        raise WPSError('Status for job "{job_id}" does not exist', job_id=job_id)

    return http.HttpResponse(job.report, content_type='text/xml')

@require_http_methods(['POST'])
@csrf_exempt
def execute(request):
    data = helpers.decoder(request.body)

    backend = backends.Backend.get_backend('CDAT')

    backend.execute(**data)

    return http.HttpResponse()
