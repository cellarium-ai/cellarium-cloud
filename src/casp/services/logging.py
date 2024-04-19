import json
import logging
import re
import typing as t
from copy import copy

import uvicorn
from starlette_context import context
from starlette_context.errors import ContextDoesNotExistError

from casp.services import settings

ANONYMOUS_USER_ID_STR: str = ""
ANONYMOUS_USER_EMAIL_STR: str = "anonymous"
UNKNOWN_USER_AGENT: str = "unknown"
UNKNOWN_TIMING: float = -1

logger = logging.getLogger(__name__)

# def setup_logging(log_level: str = "info"):
#     # remove every other logger's handlers
#     # and propagate to root logger
#     # logging.root.setLevel(log_level.upper())
#     for name in logging.root.manager.loggerDict.keys():
#         print(name)
#         # logging.getLogger(name).handlers = []
#         # logging.getLogger(name).propagate = True


class CloudTraceContext(t.NamedTuple):
    """
    Object representing the Cloud Trace context header.  The format is:
    ``
    X-Cloud-Trace-Context: TRACE_ID/SPAN_ID;o=OPTIONS
    ``
    See https://cloud.google.com/trace/docs/trace-context
    """

    trace_id: str
    span_id: str
    options: t.Optional[str]

    def get_trace(self, project_id: str) -> str:
        """
        Returns the trace field to log to Stackdriver.
        """
        return f"projects/{project_id}/traces/{self.trace_id}"

    @staticmethod
    def from_header(header: str) -> "CloudTraceContext":
        """
        Parse a Cloud Trace context header into a CloudTraceContext object.

        :param header: The value of the X-Cloud-Trace-Context header.
        :return: A CloudTraceContext object.
        """
        match = re.match(r"([a-f0-9]+)/([a-f0-9]+)(;o=(.*))?", header)
        if not match:
            logger.warning(f"Invalid X-Cloud-Trace-Context header: {header}")
        trace_id, span_id, _, options = match.groups()
        return CloudTraceContext(trace_id, span_id, options)


class StreamHandler(logging.StreamHandler):
    """
    Custom StreamHandler that logs messages as JSON objects if JSON logging is turned.
    This will also add the trace header to the log message if it is present in the context
    """

    def __init__(self, stream=None) -> None:
        super().__init__(stream=stream)
        self.log_as_json = settings.LOG_AS_JSON

    def format(self, record):
        message = super().format(record)
        if self.log_as_json:
            PROJECT = settings.GOOGLE_ACCOUNT_CREDENTIALS["project_id"]

            # Build structured log messages as an object.
            global_log_fields = {}

            # Add log correlation to nest all log messages.
            # This is only relevant in HTTP-based contexts, and is ignored elsewhere.
            # (In particular, non-HTTP-based Cloud Functions.)
            # request_is_defined = "request" in globals() or "request" in locals()
            try:
                trace_header = context.get("X-Cloud-Trace-Context")

            except ContextDoesNotExistError:
                # Logging outside of a request context.
                trace_header = None

            if trace_header and PROJECT:
                trace = CloudTraceContext.from_header(trace_header)
                global_log_fields["logging.googleapis.com/trace"] = trace.get_trace(PROJECT)
                global_log_fields["logging.googleapis.com/spanId"] = trace.span_id
                # trace = trace_header.split("/")
                # global_log_fields["logging.googleapis.com/trace"] = f"projects/{PROJECT}/traces/{trace[0]}"
                # global_log_fields["logging.googleapis.com/trace"] = f"projects/{PROJECT}/traces/{trace[0]}"

            # Complete a structured log entry.
            entry = dict(
                severity=record.levelname,
                message=message,
                # Log viewer accesses 'component' as jsonPayload.component'.
                component="arbitrary-property-for-testing",
                **global_log_fields,
            )

            return json.dumps(entry)


class CustomDefaultFormatter(uvicorn.logging.DefaultFormatter):
    """
    Custom formatter for default logs.
    """

    def __init__(self, fmt: str, fmt_json: t.Optional[str], use_colors: bool = True, *args, **kwargs):
        if settings.LOG_AS_JSON:
            # Use the JSON format if the setting is enabled.
            if fmt_json:
                fmt = fmt_json
            # Turn off colorization when logging as JSON to avoid getting control characters in the output.
            use_colors = False

        super().__init__(*args, **kwargs, use_colors=use_colors, fmt=fmt)


class CustomAccessFormatter(uvicorn.logging.AccessFormatter):
    """
    Handles the formatting of access logs and adds data to the log record if it is present in the context.  Specifically:
    - user: A loggable representation of the user from the request.
    - user_id: The internal ID of the user from the request.
    - user_email: The email address of the user from the request.
    - user_agent: The user agent string from the request.
    - request_time: The time taken to process the request. Note: the number may be slightly smaller than the actual time but should give us
                    reasonable accuracy for endpoint performance over time.
    """

    def __init__(self, fmt: str, fmt_json: t.Optional[str], use_colors: bool = True, *args, **kwargs):
        if settings.LOG_AS_JSON:
            # Use the JSON format if the setting is enabled.
            if fmt_json:
                fmt = fmt_json
            # Turn off colorization when logging as JSON to avoid getting control characters in the output.
            use_colors = False
        super().__init__(*args, **kwargs, use_colors=use_colors, fmt=fmt)

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
