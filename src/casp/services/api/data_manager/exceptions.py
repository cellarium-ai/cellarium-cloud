class DataManagerBaseException(Exception):
    """Base exception for data manager service."""


class NotFound(DataManagerBaseException):
    pass
