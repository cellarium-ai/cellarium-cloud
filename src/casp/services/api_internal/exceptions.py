class BaseAPIInternalException(Exception):
    pass


class DatabaseError(BaseAPIInternalException):
    pass


class ModelUniqueConstraintError(DatabaseError):
    pass


class IndexUniqueConstraintError(DatabaseError):
    pass


class ModelNotFoundError(DatabaseError):
    pass
