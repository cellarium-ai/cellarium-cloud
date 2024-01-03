from fastapi import Request, responses

from casp.services.api_internal.exceptions import DatabaseError


async def database_unique_constraint_error_handler(_: Request, exc: DatabaseError) -> responses.JSONResponse:
    """
    Handler for DatabaseUniqueConstraintError exceptions

    :param _:
    :param exc: DataBaseUniqueConstraintError exception

    :return: Response with 400 status code and error message
    """
    return responses.JSONResponse(
        status_code=400,
        content={"detail": str(exc)},
    )
