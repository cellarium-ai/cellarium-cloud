class APIServiceBaseException(Exception):
    """Base class for all API service exceptions"""


class AccessDeniedError(APIServiceBaseException):
    pass
