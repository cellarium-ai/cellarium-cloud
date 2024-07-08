"""
Tests things in the cellarium_general_service that can reasonably be tested with a unit test.
"""

import pytest
from mockito import mock, unstub

from casp.services import constants
from casp.services.api.services.cellarium_general_service import CellariumGeneralService
from casp.services.api.services.exceptions import InvalidClientVersionException


class TestCellariumGeneralService:

    def setup_method(self) -> None:
        self.cellarium_general_service = CellariumGeneralService(cellarium_general_dm=mock())

    def teardown_method(self) -> None:
        unstub()

    def test_validate_client_version_valid(self):
        assert self.cellarium_general_service.validate_client_version(constants.MIN_CLIENT_VERSION)

    def test_validate_client_version_too_old(self):
        assert not self.cellarium_general_service.validate_client_version("1.3.0")

    def test_validate_client_version_invalid(self):
        with pytest.raises(InvalidClientVersionException):
            self.cellarium_general_service.validate_client_version("invalid_version")
