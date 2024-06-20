class DataManagerBaseException(Exception):
    """
    Base exception for data manager service.

    These exceptions should be caught and properly marshalled into http responses.
    """


class NotFound(DataManagerBaseException):
    """
    Exception raised when a requested resource is not found in a data manager call.
    """


class CellMetadataDatabaseError(DataManagerBaseException):
    """
    Exception raised when there is an error in the cell metadata database.
    """
