import typing as t
from unittest.mock import Mock, patch

import pytest


@pytest.fixture(scope="session", autouse=True)
def patch_google_service_credentials():
    """
    Patch the `get_google_service_credentials` function globally for all tests.
    """
    mock_credentials = Mock()  # Use mock credentials
    mock_project_id = "mock-project-id"

    mock_credentials.token = "mock-token"
    mock_credentials.project_id = mock_project_id
    mock_credentials.get_project_id.return_value = mock_project_id

    def mocked_get_google_service_credentials():
        return mock_credentials, mock_project_id

    # Target the function to be patched
    patch_target = "casp.services.utils.get_google_service_credentials"

    with patch(patch_target, side_effect=mocked_get_google_service_credentials):
        yield


@pytest.fixture()
def mock_bigquery_client() -> Mock:
    return Mock()


@pytest.fixture
def patch_bigquery_client(mock_bigquery_client: Mock) -> t.Generator[None, None, None]:
    """
    Patch the BigQuery Client globally for all tests to avoid authentication.
    """

    # Patch the client class to always return the mock instance
    patch_target = "google.cloud.bigquery.Client"
    with patch(patch_target, return_value=mock_bigquery_client):
        yield
