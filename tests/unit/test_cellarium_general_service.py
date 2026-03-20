"""
Tests things in the cellarium_general_service that can reasonably be tested with a unit test.
"""

from mockito import mock, unstub
import pytest

from cellarium.cas_backend.apps.compute.services.cellarium_general_service import CellariumGeneralService
from cellarium.cas_backend.apps.compute.services.exceptions import InvalidClientVersionException
from cellarium.cas_backend.core import constants


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
