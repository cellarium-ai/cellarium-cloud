from fastapi import Request, responses, status

from casp.services._auth.exceptions import TokenException
from casp.services.api.data_manager.exceptions import NotFound
from casp.services.api.services.exceptions import APIBaseException


async def api_base_exception_handler(_: Request, exc: APIBaseException) -> responses.JSONResponse:
    """
    Handler for APIBaseException which sets the code from the exception type.

    :param _:
    :param exc: APIBaseException exception

    :return: Response with status code from the exception and its error message.
    """
    return responses.JSONResponse(
        status_code=exc.http_code,
        content={"detail": str(exc)},
    )


async def not_found_error_handler(_: Request, exc: NotFound) -> responses.JSONResponse:
    """
    Handler for NotFound exceptions

    :param _:
    :param exc: NotFound exception

    :return: Response with 400 status code and error message
    """
    return responses.JSONResponse(
        status_code=404,
        content={"detail": str(exc)},
    )


async def token_exception_handler(_: Request, exc: TokenException) -> responses.JSONResponse:
    """
    Handler for TokenException which sets the code from the exception type.

    :param _:
    :param exc: TokenException exception

    :return: Response with status code from the exception and its error message.
    """
    return responses.JSONResponse(
        status_code=exc.http_code,
        content={"detail": str(exc)},
        headers={"WWW-Authenticate": "Bearer"},
    )


async def global_error_handler(*args, **kwarg) -> responses.JSONResponse:
    """
    Top level error handler for all exceptions if not caught by more specific handlers.

    :param _:
    :param exc: Exception exception

    :return: Response with 500 status code and error message.
    """
    return responses.JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={"detail": "Error processing request"},
    )
