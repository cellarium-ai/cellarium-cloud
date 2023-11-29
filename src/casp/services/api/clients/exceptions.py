class CellariumClientBaseException(Exception):
    """Base class for all Cellarium API exceptions"""


class ClientError(CellariumClientBaseException):
    pass


class ClientTimeoutError(ClientError):
    pass


class HTTPResponseStatusCodeError(CellariumClientBaseException):
    pass


class HTTPError5XX(HTTPResponseStatusCodeError):
    pass


class HTTPError403(HTTPResponseStatusCodeError):
    pass


class HTTPError401(HTTPResponseStatusCodeError):
    pass
