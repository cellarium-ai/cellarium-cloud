
import logging
from copy import copy

from starlette_context import context
import uvicorn


ANONYMOUS_USER_ID_STR: str = ""
ANONYMOUS_USER_EMAIL_STR: str = "anonymous"
UNKNOWN_USER_AGENT: str = "unknown"
UNKNOWN_TIMING: float = -1


class CustomAccessHandler(logging.StreamHandler):
    def emit(self, record: logging.LogRecord) -> None:
        return super().emit(record)


class CustomAccessFormatter(uvicorn.logging.AccessFormatter):
    """
    Adds data to the log record if it is present in the context.  Specifically
    - user: A loggable representation of the user from the request.
    - user_id: The internal ID of the user from the request.
    - user_email: The email address of the user from the request.
    - user_agent: The user agent string from the request.
    - request_time: The time taken to process the request. Note: the number may be slightly smaller than the actual time but should give us
                    reasonable accuracy for endpoint performance over time.
    """

    def formatMessage(self, record: logging.LogRecord) -> str:
        recordcopy = copy(record)
        # Add user, user_id and user_email to the log record if they are present in the context.
        user_obj = context.data["user"] if "user" in context.data else None
        if user_obj:
            user_id = user_obj.id
            user_email = user_obj.email
            user = f"{user_id} {user_email}"
        else:
            user_id = ANONYMOUS_USER_ID_STR
            user_email = ANONYMOUS_USER_EMAIL_STR
            user = ANONYMOUS_USER_EMAIL_STR

        # Add the user agent to the log record if it is present in the context.
        if "User-Agent" in context.data:
            user_agent = context.data["User-Agent"]
        else:
            user_agent = UNKNOWN_USER_AGENT

        # Add the request time to the log record if it is present in the context.
        if "request_time" in context.data:
            request_time = context.data["request_time"]
        else:
            request_time = UNKNOWN_TIMING

        recordcopy.__dict__.update(
            {
                "user": user,
                "user_id": user_id,
                "user_email": user_email,
                "user_agent": user_agent,
                "request_time": request_time
            }
        )
        return super().formatMessage(recordcopy)
