import re

import pytest
from parameterized import parameterized

from casp.services.api.services.exceptions import InvalidInputError
from casp.services.logging import CloudTraceContext


class TestLogging:
    @parameterized.expand(
        [
            (
                "abc123/def456;o=traceoptions",
                CloudTraceContext(trace_id="abc123", span_id="def456", options="traceoptions"),
            ),
            ("abc123/def456", CloudTraceContext(trace_id="abc123", span_id="def456", options=None)),
        ]
    )
    def test_cloud_trace_context_parse(self, trace_context: str, expected: CloudTraceContext) -> None:
        assert CloudTraceContext.from_header(trace_context) == expected

    @parameterized.expand(
        [
            ("foo"),
            ("abc123nothex/def456;o=traceoptions"),
            ("abc123/def456nothex;o=traceoptions"),
            ("abc123/;o=traceoptions"),
            ("/;o=traceoptions"),
            ("/def456;o=traceoptions"),
        ]
    )
    def test_cloud_trace_context_parse_bad_trace_context(self, trace_context: str) -> None:
        with pytest.raises(InvalidInputError, match=re.escape("Invalid X-Cloud-Trace-Context header")):
            CloudTraceContext.from_header(trace_context)
