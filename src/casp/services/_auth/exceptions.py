class AuthException(Exception):
    """
    Base Auth Exception
    """

    pass


class TokenExpired(AuthException):
    pass


class TokenInvalid(AuthException):
    pass
