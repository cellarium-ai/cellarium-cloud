from fastapi import status


class APIBaseException(Exception):
    """Base class for all API service exceptions"""

    http_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR


class APIInternalError(APIBaseException):
    http_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR


class InvalidInputError(APIBaseException):
    http_code: int = status.HTTP_400_BAD_REQUEST


class AccessDeniedError(APIBaseException):
    http_code: int = status.HTTP_403_FORBIDDEN


class VectorSearchResponseError(APIInternalError):
    http_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR


class CellMetadataColumnDoesNotExist(InvalidInputError):
    http_code: int = status.HTTP_400_BAD_REQUEST


class QuotaExceededException(APIBaseException):
    http_code: int = status.HTTP_403_FORBIDDEN
