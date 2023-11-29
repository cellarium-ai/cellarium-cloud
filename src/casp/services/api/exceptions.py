class CellariumAPIBaseException(Exception):
    """Base class for all Cellarium API exceptions"""


class AnnotationServiceException(CellariumAPIBaseException):
    pass


class ServiceAPIException(CellariumAPIBaseException):
    pass
