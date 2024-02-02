class APIBaseException(Exception):
    """Base class for all API service exceptions"""


class APIInternalError(APIBaseException):
    pass


class InvalidInputError(APIBaseException):
    pass


class AccessDeniedError(APIBaseException):
    pass


class VectorSearchResponseError(APIInternalError):
    pass


class CellMetadataColumnDoesntExist(InvalidInputError):
    pass
