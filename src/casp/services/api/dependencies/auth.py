from fastapi import Depends
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sentry_sdk import set_user
from starlette_context import context

from casp.services import _auth
from casp.services.constants import ContextKeys
from casp.services.db import models

auth_scheme = HTTPBearer(auto_error=False)


async def authenticate_user(auth_token_scheme: HTTPAuthorizationCredentials = Depends(auth_scheme)) -> models.User:
    """
    Look at `Authorization` header and retrieve token by Bearer value.
    """
    if auth_token_scheme is None:
        raise _auth.exceptions.TokenInvalid("No token provided")

    jwt_token = auth_token_scheme.credentials
    user = _auth.authenticate_user_with_jwt(jwt_token)

    # Set the starlette context to know about the user making the request for things like logging
    # and quota enforcement
    context[ContextKeys.user] = user
    # Set the sentry user context
    client_host = context.get(ContextKeys.client).host if ContextKeys.client in context else None
    set_user({"id": user.id, "username": user.username, "email": user.email, "client_host": client_host})
    return user
