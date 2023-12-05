from fastapi import Request, responses

from casp.services.api.services.exceptions import AccessDeniedError, InvalidInputError


async def access_denied_error_handler(_: Request, exc: AccessDeniedError) -> responses.JSONResponse:
    """
    Handler for AccessDeniedError exceptions

    :param _:
    :param exc: AccessDeniedError exception

    :return: Response with 403 status code and error message
    """
    return responses.JSONResponse(
        status_code=403,
        content={"detail": str(exc)},
    )


async def invalid_input_error_handler(_: Request, exc: InvalidInputError) -> responses.JSONResponse:
    """
    Handler for InvalidInputError exceptions

    :param _:
    :param exc: InvalidInputError exception

    :return: Response with 400 status code and error message
    """
    return responses.JSONResponse(
        status_code=400,
        content={"detail": str(exc)},
    )
