from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sentry_sdk import set_user
from starlette_context import context

from casp.services import _auth
from casp.services.constants import ContextKeys
from casp.services.db import models

auth_scheme = HTTPBearer()


async def authenticate_user(auth_token_scheme: HTTPAuthorizationCredentials = Depends(auth_scheme)) -> models.User:
    """
    Look at `Authorization` header and retrieve token by Bearer value.
    """
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    jwt_token = auth_token_scheme.credentials
    try:
        user = _auth.authenticate_user_with_jwt(jwt_token)
    except _auth.exceptions.TokenInvalid or _auth.exceptions.TokenExpired:
        raise credentials_exception

    # Set the starlette context to know about the user making the request for things like logging
    # and quota enforcement
    context[ContextKeys.user] = user
    # Set the sentry user context
    client_host = context.get(ContextKeys.client).host if ContextKeys.client in context else None
    set_user({"id": user.id, "username": user.username, "email": user.email, "client_host": client_host})
    return user
