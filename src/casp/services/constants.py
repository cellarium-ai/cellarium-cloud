from enum import Enum


class ContextKeys(str, Enum):
    """
    Keys used in the context object to store information about the request.
    """

    user = "user"
    client = "client"


class HeaderKeys(str, Enum):
    """
    Keys used in the headers of the request.
    """

    cloud_trace_context = "X-Cloud-Trace-Context"
    authorization = "Authorization"


class LogRecordKeys(str, Enum):
    """
    Keys used in the log records.
    """

    user = ("user",)
    user_id = ("user_id",)
    user_email = ("user_email",)
    user_agent = ("user_agent",)
    request_time = ("request_time",)
