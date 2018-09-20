import datetime
import collections
import json
import logging
import re
import random
import string

from django import http
from django import db
from django.conf import settings
from django.core.mail import send_mail
from django.contrib.auth import authenticate
from django.contrib.auth import login
from django.contrib.auth import logout
from django.shortcuts import redirect
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.decorators.http import require_http_methods
from django.views.decorators.cache import never_cache
from openid.consumer import discover
from myproxy.client import MyProxyClient

from wps import forms
from wps import metrics
from wps import models
from wps import WPSError
from wps.auth import oauth2
from wps.auth import openid
from wps.views import common

logger = logging.getLogger('wps.views.auth')

URN_AUTHORIZE = 'urn:esg:security:oauth:endpoint:authorize'
URN_ACCESS = 'urn:esg:security:oauth:endpoint:access'
URN_RESOURCE = 'urn:esg:security:oauth:endpoint:resource'
URN_MPC = 'urn:esg:security:myproxy-service'

discover.OpenIDServiceEndpoint.openid_type_uris.extend([
    URN_AUTHORIZE,
    URN_ACCESS,
    URN_RESOURCE,
    URN_MPC
])

FORGOT_USERNAME_SUBJECT = 'CWT WPS Username Recovery'
FORGOT_USERNAME_MESSAGE = """
Hello {username},

You've requested the recovery of your username: {username}.

Thank you,
ESGF CWT Team
"""

FORGOT_PASSWORD_SUBJECT = 'CWT WPS Password Reset'
FORGOT_PASSWORD_MESSAGE = """
Hello {username},
<br><br>
You've request the reset of you password. Please follow this <a href="{reset_url}">link</a> to reset you password.
<br><br>
Thank you,
ESGF CWT Team
"""

CREATE_SUBJECT = 'Welcome to ESGF compute server'
CREATE_MESSAGE = """
Thank you for creating an account for the ESGF compute server. Please login into your account <a href="{login_url}">here</a>.

If you have any questions or concerns please email the <a href="mailto:{admin_email}">server admin</a>.
"""

class ResetPasswordInvalidStateError(WPSError):
    def __init__(self):
        msg = 'Invalid state while recovering password, please try again'

        super(ResetPasswordInvalidStateError, self).__init__(msg)

class ResetPasswordTokenExpiredError(WPSError):
    def __init__(self):
        msg = 'Token to reset password has expired'

        super(ResetPasswordTokenExpiredError, self).__init__(msg)

class ResetPasswordTokenMismatchError(WPSError):
    def __init__(self):
        msg = 'Token to reset password does not match expected value'

        super(ResetPasswordTokenMismatchError, self).__init__(msg)

class MPCEndpointParseError(WPSError):
    def __init__(self):
        msg = 'Parsing host/port from OpenID services failed'

        super(MPCEndpointParseError, self).__init__(msg)

@require_http_methods(['GET'])
@ensure_csrf_cookie
def user_cert(request):
    try:
        if not settings.CERT_DOWNLOAD_ENABLED:
            return http.HttpResponseBadRequest()

        common.authentication_required(request)

        user = request.user

        cert = user.auth.cert

        content_type = 'application/force-download' 

        response = http.HttpResponse(cert, content_type=content_type)

        response['Content-Disposition'] = 'attachment; filename="cert.pem"'

        response['Content-Length'] = len(cert)
    except WPSError as e:
        return http.HttpResponseBadRequest()
    else:
        return response

@require_http_methods(['POST'])
@ensure_csrf_cookie
def create(request):
    try:
        form = forms.CreateForm(request.POST)

        data = common.validate_form(form, ('username', 'email', 'openid', 'password'))

        try:
            user = models.User.objects.create_user(data['username'], data['email'], data['password'])
        except db.IntegrityError:
            raise common.DuplicateUserError(username=data['username'])

        models.Auth.objects.create(openid_url=data['openid'], user=user)

        try:
            send_mail(CREATE_SUBJECT,
                      '',
                      settings.WPS_ADMIN_EMAIL,
                      [user.email],
                      html_message=CREATE_MESSAGE.format(login_url=settings.WPS_LOGIN_URL, admin_email=settings.WPS_ADMIN_EMAIL))
        except Exception:
            logger.exception('Error sending confirmation email')
            
            pass
    except WPSError as e:
        logger.exception('Error creating account')

        return common.failed(str(e))
    else:
        return common.success('Successfully created account for "{}"'.format(data['username']))

@require_http_methods(['POST'])
@ensure_csrf_cookie
def user_login_openid(request):
    try:
        form = forms.OpenIDForm(request.POST)

        data = common.validate_form(form, ('openid_url',))

        url = openid.begin(request, data['openid_url'])
    except WPSError as e:
        logger.exception('Error logging user in with OpenID')

        return common.failed(str(e))
    else:
        return common.success({'redirect': url})

@require_http_methods(['GET'])
@ensure_csrf_cookie
def user_login_openid_callback(request):
    try:
        openid_url, attrs = openid.complete(request)

        try:
            user = models.User.objects.get(auth__openid_url=openid_url)
        except models.User.DoesNotExist:
            username = openid_url.split('/')[-1]

            user = models.User.objects.create_user(username, attrs['email'])

            models.Auth.objects.create(openid_url=openid_url, user=user)

        login(request, user)
    except WPSError as e:
        logger.exception('Error handling OpenID callback')

        return common.failed(str(e))
    else:
        return redirect('{}?expires={}'.format(settings.WPS_OPENID_CALLBACK_SUCCESS, request.session.get_expiry_date()))

@require_http_methods(['POST'])
@ensure_csrf_cookie
def user_login(request):
    try:
        form = forms.LoginForm(request.POST)

        data = common.validate_form(form, ('username', 'password'))

        logger.info('Attempting to login user {}'.format(data['username']))

        user = authenticate(request, username=data['username'], password=data['password'])

        if user is None:
            raise common.AuthenticationError(user=data['username'])

        login(request, user)
    except WPSError as e:
        logger.exception('Error logging user in')

        return common.failed(str(e))
    else:
        data = common.user_to_json(user)

        data['expires'] = request.session.get_expiry_date()

        return common.success(data)

@require_http_methods(['GET'])
@ensure_csrf_cookie
def user_logout(request):
    try:
        common.authentication_required(request)

        logger.info('Logging user {} out'.format(request.user.username))

        logout(request)
    except WPSError as e:
        logger.exception('Error logging user out')

        return common.failed(str(e))
    else:
        return common.success('Logged out')

@require_http_methods(['POST'])
@ensure_csrf_cookie
def login_oauth2(request):
    try:
        common.authentication_required(request)

        logger.info('Authenticating OAuth2 for {}'.format(request.user.auth.openid_url))

        auth_service, cert_service = openid.services(request.user.auth.openid_url, (URN_AUTHORIZE, URN_RESOURCE))

        redirect_url, state = oauth2.get_authorization_url(auth_service.server_url, cert_service.server_url)

        logger.info('Retrieved authorization url for OpenID {}'.format(request.user.auth.openid_url))

        request.session.update({
            'oauth_state': state,
            'openid': request.user.auth.openid_url
        })
    except WPSError as e:
        logger.exception('Error authenticating OAuth2')

        return common.failed(str(e))
    else:
        return common.success({'redirect': redirect_url})

@require_http_methods(['GET'])
@ensure_csrf_cookie
def oauth2_callback(request):
    user = None

    try:
        openid_url = request.session.pop('openid')

        oauth_state = request.session.pop('oauth_state')

        logger.info('Handling OAuth2 callback for %r current state %r',
                    openid_url, oauth_state)

        user = models.User.objects.get(auth__openid_url = openid_url)

        logger.info('Discovering token and certificate services')

        token_service, cert_service = openid.services(openid_url, (URN_ACCESS, URN_RESOURCE))

        request_url = '{}?{}'.format(settings.WPS_OAUTH2_CALLBACK, request.META['QUERY_STRING'])

        logger.info('Getting token from service')

        token = oauth2.get_token(token_service.server_url, request_url, oauth_state)

        logger.info('Getting certificate from service')

        cert, key, new_token = oauth2.get_certificate(token, oauth_state, token_service.server_url, cert_service.server_url)

        logger.info('Updating user with token, certificate and state')

        user.auth.update('oauth2', [cert, key], token=new_token, state=oauth_state)
    except KeyError as e:
        logger.exception('Missing %r key from session data', e)

        pass
    except (WPSError, oauth2.OAuth2Error) as e:
        logger.exception('OAuth2 callback failed')

        if user is not None:
            extra = json.loads(user.auth.extra)

            extra['error'] = 'OAuth2 callback failed "{}"'.format(str(e))

            user.auth.extra = json.dumps(extra)

            user.auth.save()
    
    logger.info('Finished handling OAuth2 callback, redirect to profile')

    return redirect(settings.WPS_PROFILE_URL)

@require_http_methods(['POST'])
@ensure_csrf_cookie
def login_mpc(request):
    try:
        common.authentication_required(request)

        form = forms.MPCForm(request.POST)

        data = common.validate_form(form, ('username', 'password'))

        logger.info('Authenticating MyProxyClient for {}'.format(data['username']))

        services = openid.services(request.user.auth.openid_url, (URN_MPC,))

        g = re.match('socket://(.*):(.*)', services[0].server_url)

        if g is None or len(g.groups()) != 2:
            raise MPCEndpointParseError()

        host, port = g.groups()

        from OpenSSL import SSL

        MyProxyClient.SSL_METHOD = SSL.TLSv1_2_METHOD

        try:
            m = MyProxyClient(hostname=host, caCertDir=settings.WPS_CA_PATH)

            c = m.logon(data['username'], data['password'], bootstrap=True)
        except Exception as e:
            raise common.AuthenticationError(user=data['username'])

        logger.info('Authenticated with MyProxyClient backend for user {}'.format(data['username']))

        request.user.auth.update('myproxyclient', c)
    except WPSError as e:
        logger.exception('Error authenticating MyProxyClient')

        return common.failed(str(e))
    else:
        return common.success({
            'type': request.user.auth.type,
            'api_key': request.user.auth.api_key
        })

@require_http_methods(['GET'])
@ensure_csrf_cookie
def forgot_username(request):
    try:
        try:
            email = request.GET['email']
        except KeyError as e:
            raise common.MissingParameterError(name=str(e))

        logger.info('Recovering username for "{}"'.format(email))

        try:
            user = models.User.objects.get(email=email)
        except models.User.DoesNotExist:
            raise common.UserEmailDoesNotExistError(email=email)

        try:
            send_mail(FORGOT_USERNAME_SUBJECT,
                      FORGOT_USERNAME_MESSAGE.format(username=user.username, login_url=settings.WPS_LOGIN_URL),
                      settings.WPS_ADMIN_EMAIL,
                      [user.email])
        except:
            raise common.MailError(email=user.email)
    except WPSError as e:
        logger.exception('Error recovering username')

        return common.failed(str(e))
    else:
        return common.success({'redirect': settings.WPS_LOGIN_URL})

@require_http_methods(['GET'])
@ensure_csrf_cookie
def forgot_password(request):
    try:
        try:
            username = request.GET['username']
        except KeyError as e:
            raise common.MissingParameterError(name=str(e))

        logger.info('Starting password reset process for user "{}"'.format(username))

        try:
            user = models.User.objects.get(username=username)
        except models.User.DoesNotExist:
            raise common.UserDoesNotExistError(username=username)

        try:
            extra = json.loads(user.auth.extra)
        except Exception:
            extra = {}

        extra['reset_token'] = ''.join(random.choice(string.ascii_letters + string.digits) for _ in xrange(64))

        extra['reset_expire'] = datetime.datetime.now() + datetime.timedelta(1)

        user.auth.extra = json.dumps(extra, default=lambda x: x.strftime('%x %X'))

        user.auth.save()

        reset_url = '{}?token={}&username={}'.format(settings.WPS_PASSWORD_RESET_URL, extra['reset_token'], user.username)

        try:
            send_mail(FORGOT_PASSWORD_SUBJECT,
                      '',
                      settings.WPS_ADMIN_EMAIL,
                      [user.email],
                      html_message=FORGOT_PASSWORD_MESSAGE.format(username=user.username,
                                                              reset_url=reset_url)
                      )
        except:
            raise common.MailError(email=user.email)
    except WPSError as e:
        return common.failed(str(e))
    else:
        return common.success({'redirect': settings.WPS_LOGIN_URL})

@require_http_methods(['GET'])
@ensure_csrf_cookie
def reset_password(request):
    try:
        try:
            token = str(request.GET['token'])

            username = str(request.GET['username'])

            password = str(request.GET['password'])
        except KeyError as e:
            raise common.MissingParameterError(name=str(e))

        logger.info('Resetting password for "{}"'.format(username))

        try:
            user = models.User.objects.get(username=username)
        except models.User.DoesNotExist:
            raise common.UserDoesNotExistError(username=username)

        try:
            extra = json.loads(user.auth.extra)
        except Exception:
            extra = {}

        try:
            reset_token = extra['reset_token']

            reset_expire = extra['reset_expire']
        except KeyError as e:
            raise ResetPasswordInvalidStateError()
        else:
            del extra['reset_token']

            del extra['reset_expire']

        expires = datetime.datetime.strptime(reset_expire, '%x %X')

        if datetime.datetime.now() > expires:
            raise ResetPasswordTokenExpiredError()

        if reset_token != token:
            raise ResetPasswordTokenMismatchError()

        user.auth.extra = json.dumps(extra)

        user.auth.save()

        user.set_password(password)

        user.save()

        logger.info('Successfully reset password for "{}"'.format(username))
    except WPSError as e:
        return common.failed(str(e))
    else:
        return common.success({'redirect': settings.WPS_LOGIN_URL})
