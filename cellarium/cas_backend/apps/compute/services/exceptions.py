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
    http_code: int = status.HTTP_429_TOO_MANY_REQUESTS


class InvalidClientVersionException(APIBaseException):
    http_code: int = status.HTTP_400_BAD_REQUEST

    def __init__(self, client_version: str):
        super().__init__(f"Client version {client_version} is not a valid version string")


class InvalidAnndataFileException(InvalidInputError):
    http_code: int = status.HTTP_400_BAD_REQUEST
