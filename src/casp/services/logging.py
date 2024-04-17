import json
import logging
from copy import copy

import uvicorn
from starlette_context import context
from starlette_context.errors import ContextDoesNotExistError

from casp.services import settings

ANONYMOUS_USER_ID_STR: str = ""
ANONYMOUS_USER_EMAIL_STR: str = "anonymous"
UNKNOWN_USER_AGENT: str = "unknown"
UNKNOWN_TIMING: float = -1


class InterceptHandler(logging.Handler):
    def emit(self, record):
        # Get corresponding Loguru level if it exists
        # try:
        #     level = logger.level(record.levelname).name
        # except ValueError:
        # level = record.levelno

        # Find caller from where originated the logged message
        # frame, depth = logging.currentframe(), 2
        # while frame.f_code.co_filename == logging.__file__:
        #     frame = frame.f_back
        #     depth += 1

        # logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())
        # Uncomment and populate this variable in your code:
        PROJECT = settings.GOOGLE_ACCOUNT_CREDENTIALS["project_id"]

        # Build structured log messages as an object.
        global_log_fields = {}

        # Add log correlation to nest all log messages.
        # This is only relevant in HTTP-based contexts, and is ignored elsewhere.
        # (In particular, non-HTTP-based Cloud Functions.)
        # request_is_defined = "request" in globals() or "request" in locals()
        # if request_is_defined and request:
        #     trace_header = request.headers.get("X-Cloud-Trace-Context")
        try:
            trace_header = context.get("X-Cloud-Trace-Context")
        except ContextDoesNotExistError:
            # Logging outside of a request context.
            trace_header = None

        # if "X-Cloud-Trace-Context" in context.data else None

        if trace_header and PROJECT:
            trace = trace_header.split("/")
            global_log_fields["logging.googleapis.com/trace"] = f"projects/{PROJECT}/traces/{trace[0]}"

        # Complete a structured log entry.
        entry = dict(
            severity=record.levelname,
            message=record.getMessage(),
            # Log viewer accesses 'component' as jsonPayload.component'.
            # component="arbitrary-property",
            **global_log_fields,
        )

        print(json.dumps(entry))


def setup_logging(log_level: str = "info", json_logs: bool = True):
    if json_logs:
        # intercept everything at the root logger
        logging.root.handlers = [InterceptHandler()]
        logging.root.setLevel(log_level.upper())

        # remove every other logger's handlers
        # and propagate to root logger
        for name in logging.root.manager.loggerDict.keys():
            logging.getLogger(name).handlers = []
            logging.getLogger(name).propagate = True


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
                "request_time": request_time,
            }
        )
        return super().formatMessage(recordcopy)


class JsonHandler(logging.StreamHandler):
    def emit(self, record: logging.LogRecord) -> None:
        return super().emit(record)
