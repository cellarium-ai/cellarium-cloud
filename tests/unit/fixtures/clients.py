"""
Fixtures for mocking the MatchingClient used in unit tests.

This module provides pytest fixtures for mocking the behavior of the `MatchingClient`. These fixtures
allow unit tests to simulate interactions with the MatchingClient without requiring an actual backend.
"""

import typing as t
from unittest.mock import patch

import pytest

from casp.services.api import schemas
from tests.unit.fixtures import constants, mocks


@pytest.fixture
def mock_matching_client(cell_info_data: t.List[schemas.CellariumCellMetadata]) -> mocks.MockMatchingClient:
    """
    Fixture for a simplified mocked `MatchingClient` that uses the mock_cell_source.

    :param cell_info_data: Source of cell data with unique IDs and cell types.

    :return: A `MockMatchingClient` instance that generates mock neighbors.
    """
    return mocks.MockMatchingClient(cell_info_data=cell_info_data)


@pytest.fixture
def patch_matching_client(mock_matching_client) -> t.Generator[None, None, None]:
    """
    Fixture to automatically patch the `MatchingClient.from_index` method with a mocked instance.

    This fixture ensures that any test using the `MatchingClient.from_index` method will automatically
    receive the mocked instance provided by the `mock_matching_client` fixture. When this fixture is used by a test
    it will be applied before its execution.

    **Scope**: Function

    **Autouse**: No

    :param mock_matching_client: The mocked instance of `MatchingClient` provided by the `mock_matching_client` fixture.
    """

    with patch(target=constants.MATCHING_CLIENT_FROM_INDEX_TARGET, return_value=mock_matching_client):
        yield
