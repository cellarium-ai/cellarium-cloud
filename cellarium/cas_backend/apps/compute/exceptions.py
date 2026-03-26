from fastapi import status


class APIBaseException(Exception):
    """Base class for all API service exceptions"""

    http_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR


class APIInternalError(APIBaseException):
    http_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR
