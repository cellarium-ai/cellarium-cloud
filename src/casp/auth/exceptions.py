class AuthException(Exception):
    pass


class TokenExpired(AuthException):
    pass


class TokenInvalid(AuthException):
    pass
