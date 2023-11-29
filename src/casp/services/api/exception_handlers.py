from fastapi import Request, responses

from casp.services.api.main import application
from casp.services.api.services.exceptions import AccessDeniedError


@application.exception_handler(AccessDeniedError)
async def access_denied_error_handler(_: Request, exc: AccessDeniedError) -> responses.JSONResponse:
    """
    Handler for AccessDeniedError exceptions

    :param _:
    :param exc: AccessDeniedError exception

    :return: Response with 403 status code and error message
    """
    return responses.JSONResponse(
        status_code=403,
        content={"detail": exc},
    )
