from enum import Enum

from fastapi import Depends
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sentry_sdk import set_user
from starlette_context import context

from casp.services import _auth
from casp.services.constants import ContextKeys
from casp.services.db import models

auth_scheme = HTTPBearer(auto_error=False)


class TokenTypes(str, Enum):
    """
    Types of authentication tokens that can be used.
    """

    JWT = "jwt"
    OPAQUE = "opaque"

    @staticmethod
    def detect_type(token: str) -> "TokenTypes":
        """
        Detect the type of token based on the token string.

        :param token: Token string to detect the type of
        :raises TokenInvalid: If the token type is not recognized
        :return: TokenTypes enum value
        """

        if token.startswith("Bearer "):
            token = token[7:]

        # Naive detection of ot the type of auth token based on part count
        num_parts = len(token.split("."))
        if num_parts == 2:
            return TokenTypes.OPAQUE
        elif num_parts == 3:
            return TokenTypes.JWT
        else:
            raise _auth.exceptions.TokenInvalid("Token type not recognized")


async def authenticate_user(auth_token_scheme: HTTPAuthorizationCredentials = Depends(auth_scheme)) -> models.User:
    """
    Look at `Authorization` header and retrieve token by Bearer value.
    """
    if auth_token_scheme is None:
        raise _auth.exceptions.TokenInvalid("No token provided")

    token = auth_token_scheme.credentials
    token_type = TokenTypes.detect_type(token)
    match token_type:
        case TokenTypes.OPAQUE:
            user = _auth.authenticate_user_with_opaque_token(token)
        case TokenTypes.JWT:
            user = _auth.authenticate_user_with_jwt(token)
        case _:
            # This should never reach here but is left as protection in case a new token type is added
            raise _auth.exceptions.TokenInvalid("Token type not recognized")

    # Set the starlette context to know about the user making the request for things like logging
    # and quota enforcement
    context[ContextKeys.user] = user
    # Set the sentry user context
    client_host = context.get(ContextKeys.client).host if ContextKeys.client in context else None
    set_user({"id": user.id, "username": user.username, "email": user.email, "client_host": client_host})
    return user
