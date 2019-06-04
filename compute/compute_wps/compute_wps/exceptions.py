from compute_wps.util import wps_response


class WPSError(Exception):
    def __init__(self, text, *args, **kwargs):
        self.code = kwargs.get('code', None)

        if self.code is None:
            self.code = wps_response.NoApplicableCode

        super(WPSError, self).__init__(text.format(*args, **kwargs))


class AccessError(WPSError):
    def __init__(self, url, error):
        msg = 'Error accessing {!s}: {!s}'

        super(AccessError, self).__init__(msg, url, error)
