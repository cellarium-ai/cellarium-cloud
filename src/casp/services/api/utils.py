from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from casp.services import _auth
from casp.services.db import models

auth_scheme = HTTPBearer()


async def get_current_user(auth_token_scheme: HTTPAuthorizationCredentials = Depends(auth_scheme)) -> models.User:
    """
    Look at `Authorization` header and retrieve token by Bearer value.

    Could be used as a request middleware.
    """
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"Authenticate": "Bearer"},
    )
    jwt_token = auth_token_scheme.credentials
    try:
        user = _auth.authenticate_user_with_jwt(jwt_token)
    except _auth.exceptions.TokenInvalid or _auth.exceptions.TokenExpired:
        raise credentials_exception
    return user
