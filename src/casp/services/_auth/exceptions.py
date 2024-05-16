from fastapi import status


class TokenException(Exception):
    default_message: str = "Token Error"
    """
    Base Token-based Exception
    """

    def __init__(self, *args, **kwargs):
        if not args:
            args = (self.default_message,)

        super().__init__(*args, **kwargs)
        self.http_code = status.HTTP_401_UNAUTHORIZED


class TokenExpired(TokenException):
    default_message: str = "Token is Expired"


class TokenInvalid(TokenException):
    default_message: str = "Token is Invalid"


class TokenInactive(TokenException):
    default_message: str = "Token is Invalid"


class TokenGenerationError(TokenException):
    default_message: str = "Could not generate token"
