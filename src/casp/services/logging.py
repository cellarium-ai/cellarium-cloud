import json
import logging
import re
import typing as t
from copy import copy

import uvicorn
from starlette_context import context
from starlette_context.errors import ContextDoesNotExistError
from starlette_context.header_keys import HeaderKeys as BaseHeaderKeys

from casp.services import settings
from casp.services.api.services.exceptions import InvalidInputError
from casp.services.constants import ContextKeys, HeaderKeys, LogRecordKeys
from casp.services.utils import get_google_service_credentials

ANONYMOUS_USER_ID_STR: str = ""
ANONYMOUS_USER_EMAIL_STR: str = "anonymous"
UNKNOWN_USER_AGENT: str = "unknown"
UNKNOWN_TIMING: float = -1

logger = logging.getLogger(__name__)


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
        Parse a Cloud Trace context header into a CloudTraceContext object.  Throws an exception if the header is invalid.

        :param header: The value of the X-Cloud-Trace-Context header.
        :return: A CloudTraceContext object.
        """
        match = re.match(r"^([a-f0-9]+)/([a-f0-9]+)((;o=(.*))?)$", header)
        if not match:
            logger.warning(f"Invalid X-Cloud-Trace-Context header: {header}")
            raise InvalidInputError("Invalid X-Cloud-Trace-Context header")
        trace_id, span_id, _, _, options = match.groups()
        return CloudTraceContext(trace_id, span_id, options)


class StreamHandler(logging.StreamHandler):
    """
    Custom StreamHandler that logs messages as JSON objects if JSON logging is turned on.
    This will also add the trace header to the log message if it is present in the context
    """

    def __init__(self, stream=None) -> None:
        super().__init__(stream=stream)
        self.log_as_json = settings.LOG_AS_JSON
        _, self.project_id = get_google_service_credentials()

    def format(self, record):
        message = super().format(record)
        if self.log_as_json:

            # Build structured log messages as an object.
            global_log_fields = {}

            # Add log correlation to nest all log messages.
            # This is only relevant in HTTP-based contexts, and is ignored elsewhere.
            sentry_trace_id = None
            try:
                google_trace_header = context.get(HeaderKeys.cloud_trace_context)
                if ContextKeys.sentry_trace_id in context.data:
                    sentry_trace_id = context.data[ContextKeys.sentry_trace_id]

            except ContextDoesNotExistError:
                # Logging outside of a request context.
                google_trace_header = None

            if google_trace_header and self.project_id:
                trace = CloudTraceContext.from_header(google_trace_header)
                global_log_fields["logging.googleapis.com/trace"] = trace.get_trace(self.project_id)
                global_log_fields["logging.googleapis.com/spanId"] = trace.span_id

            # Complete a structured log entry.
            entry = dict(
                severity=record.levelname,
                message=message,
                # Add ability to filter logs on sentry trace id
                sentry_trace_id=sentry_trace_id,
                **global_log_fields,
            )

            return json.dumps(entry)
        else:
            return super().format(record)


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
        user_obj = context.data[ContextKeys] if ContextKeys in context.data else None
        if user_obj:
            user_id = user_obj.id
            user_email = user_obj.email
            user = f"{user_id} {user_email}"
        else:
            user_id = ANONYMOUS_USER_ID_STR
            user_email = ANONYMOUS_USER_EMAIL_STR
            user = ANONYMOUS_USER_EMAIL_STR

        # Add the user agent to the log record if it is present in the context.
        if BaseHeaderKeys.user_agent in context.data:
            user_agent = context.data[BaseHeaderKeys.user_agent]
        else:
            user_agent = UNKNOWN_USER_AGENT

        recordcopy.__dict__.update(
            {
                LogRecordKeys.user.value: user,
                LogRecordKeys.user_id.value: user_id,
                LogRecordKeys.user_email.value: user_email,
                LogRecordKeys.user_agent.value: user_agent,
            }
        )
        return super().formatMessage(recordcopy)
