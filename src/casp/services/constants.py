from enum import Enum


class ContextKeys(str, Enum):
    """
    Keys used in the context object to store information about the request.
    """

    # The user object associated with the request.
    user = "user"
    # The client object associated with the request.
    client = "client"
    # The sentry trace id associated with the request.
    sentry_trace_id = "sentry_trace_id"


class HeaderKeys(str, Enum):
    """
    Keys used in the headers of the request.
    """

    # The GCP-ingected trace context header.
    cloud_trace_context = "X-Cloud-Trace-Context"
    # The authorization header.
    authorization = "Authorization"
    # The trace id to return to users as a header.  Can be used for debugging.
    trace_id = "x-trace-id"


class LogRecordKeys(str, Enum):
    """
    Keys used in the log records.
    """

    # A representation of the the user as <id> <email> if authorized or anonymous if not or if no auth was provided
    user = "user"
    # The id of the user.
    user_id = "user_id"
    # The email of the user.
    user_email = "user_email"
    # The browser/client of the requestor.
    user_agent = "user_agent"
