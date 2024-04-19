from enum import Enum


class HeaderKeys(str, Enum):
    cloud_trace_context = "X-Cloud-Trace-Context"
    authorization = "Authorization"
