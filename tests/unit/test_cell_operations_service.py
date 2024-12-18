"""
Unit Tests for Cell Operations Service and Related Functionalities

This module provides tests for the `CellOperationsService` class, its interaction with related services, and
its handling of database models and external dependencies like the MatchingClient. It also includes
pytest fixtures to set up the necessary environment for these tests.

"""

import io
import math
import typing as t
from unittest.mock import Mock, patch

import anndata
import numpy as np
import pytest
from sqlalchemy.orm import Session
from starlette_context.ctx import _request_scope_context_storage

from casp.services import settings
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.services import exceptions
from casp.services.api.services.cell_operations_service import CellOperationsService
from casp.services.api.services.consensus_engine.strategies.ontology_aware import CellOntologyResource
from casp.services.constants import ContextKeys
from casp.services.db import models
from tests.unit.fixtures import constants, mocks


@pytest.fixture
def patch_starlette_context() -> t.Generator[t.Dict[str, str], None, None]:
    """
    Fixture to patch the `starlette_context`'s `_request_scope_context_storage` globally during tests.

    Temporarily sets mock data in the Starlette request context to avoid `LookupError`
    and simulate a test environment with predefined context values.

    **Scope**: Function

    :return: A dictionary containing the patched context data.
    """
    mock_data = {ContextKeys.sentry_trace_id: constants.TEST_REQUEST_ID_NEW}

    # Directly set the value in the `ContextVar`
    token = _request_scope_context_storage.set(mock_data)

    try:
        yield
    finally:
        # Reset the `ContextVar` after the test
        _request_scope_context_storage.reset(token)


@pytest.fixture
def cas_model(db_session: Session) -> models.CASModel:
    """
    Provides a CASModel instance from the test database.

    :param db_session: Database session fixture.
    :return: CASModel instance.
    """
    return db_session.query(models.CASModel).filter(models.CASModel.model_name == constants.TEST_MODEL_NAME).first()


@pytest.fixture
def user_with_quota(db_session: Session) -> models.User:
    """
    Fixture to provide a user with sufficient quota from the database.

    :param db_session: Fixture providing the database session.
    :return: User instance with sufficient quota.
    """
    return db_session.query(models.User).filter(models.User.username == constants.USER_EMAIL_WITH_QUOTA).first()


@pytest.fixture
def user_without_quota(db_session: Session) -> models.User:
    """
    Fixture to provide a user without sufficient quota from the database.

    :param db_session: Fixture providing the database session.
    :return: User instance without sufficient quota.
    """
    return db_session.query(models.User).filter(models.User.username == constants.USER_EMAIL_WITHOUT_QUOTA).first()


@pytest.fixture(
    params=[
        (np.random.rand(10, 32), 3, 4),  # 10 embeddings, chunk size 3
        (np.random.rand(15, 32), 5, 3),  # 15 embeddings, chunk size 5
        (np.random.rand(10, 32), 10, 1),  # Exactly one chunk
        (np.random.rand(10, 32), 20, 1),  # Chunk size larger than data
        (np.random.rand(0, 32), 5, 0),  # Empty embeddings
        (np.random.rand(11, 32), 4, 3),  # Irregular chunks
    ]
)
def embedding_test_data(request: pytest.FixtureRequest) -> t.Tuple[np.ndarray, int, int]:
    """
    Fixture to provide test data for chunk splitting tests.

    :param request: Built-in pytest fixture to iterate through parameters.
    :return: Tuple of (embeddings, chunk_size, expected_num_chunks).
    """
    return request.param


def test__split_embeddings_into_chunks(embedding_test_data: t.Tuple[np.ndarray, int, int]):
    """
    Test the static method `_split_embeddings_into_chunks` for correctness.

    Validate that the method correctly splits embeddings into the specified chunk size,
    preserving data integrity and ensuring correct chunk counts.

    :param embedding_test_data: Embedding data for the test.
    """
    embeddings, chunk_size, expected_num_chunks = embedding_test_data
    chunks = CellOperationsService._CellOperationsService__split_embeddings_into_chunks(embeddings, chunk_size)

    # Assert
    assert len(chunks) == expected_num_chunks

    # Validate the size of each chunk (all except for the last). Don't check if expected_num_chunks == 0
    if expected_num_chunks > 0:
        for chunk in chunks[:-1]:
            assert len(chunk) == chunk_size

        # Validate the size of the last chunk
        assert len(chunks[-1] <= chunk_size)

    # Validate that the data is correctly preserved across chunks
    if len(embeddings) > 0:
        reconstructed = np.vstack(chunks)
        np.testing.assert_array_equal(reconstructed, embeddings)


@pytest.mark.asyncio
async def test__get_knn_matches_for_chunk(
    populate_db: None,
    patch_bigquery_client: None,
    mock_matching_client: mocks.MockMatchingClient,
    cell_operations_service_with_mocks: CellOperationsService,
    embedding_test_data: t.Tuple[np.ndarray, int, int],
):
    """
    Test the private `_get_knn_matches_for_chunk` method.

    Validate that the method retrieves nearest neighbors for a given chunk of embeddings using the mocked
    `MatchingClient`.

    :param populate_db: Use a fixture to populate the test database.
    :param mock_matching_client: Mocked MatchingClient instance.
    :param patch_bigquery_client: Fixture to patch BigQuery client.
    :param cell_operations_service_with_mocks: Mocked CellOperationsService instance.
    :param embedding_test_data: Embedding test data.
    """
    embeddings_chunk, _, _ = embedding_test_data
    service = cell_operations_service_with_mocks

    result = await service._CellOperationsService__get_knn_matches_for_chunk(
        embeddings_chunk=embeddings_chunk, client=mock_matching_client
    )

    assert isinstance(result, MatchResult)
    assert len(result.matches) == len(embeddings_chunk)

    # Ensure `match` was called with the correct arguments
    mock_matching_client.match.assert_awaited_once_with(queries=embeddings_chunk)


@pytest.mark.asyncio
async def test_get_knn_matches_from_embeddings(
    populate_db: None,
    patch_matching_client: None,
    patch_bigquery_client: None,
    mock_matching_client: mocks.MockMatchingClient,
    cell_operations_service_with_mocks: CellOperationsService,
    embedding_test_data: t.Tuple[np.ndarray, int, int],
    cas_model: models.CASModel,
):
    """
    Test the `get_knn_matches_from_embeddings` method.

    Validates the retrieval of K-Nearest Neighbors for a batch of embeddings, ensuring correct interactions with the
    mocked `MatchingClient`.

    :param populate_db: Use a fixture to populate the test database.
    :param patch_matching_client: Use a fixture to patch the MatchingClient.
    :param patch_bigquery_client: Use a fixture to patch BigQuery client.
    :param mock_matching_client: Mocked MatchingClient instance.
    :param cell_operations_service_with_mocks: Mocked CellOperationsService instance.
    :param cas_model: CASModel instance from the test database.
    """
    # embeddings = np.random.rand(10, 64)  # A chunk of 5 embeddings, each of 64 dimensions
    embeddings, _, _ = embedding_test_data
    service = cell_operations_service_with_mocks
    result = await service.get_knn_matches_from_embeddings(embeddings=embeddings, model=cas_model)

    assert isinstance(result, MatchResult)
    assert len(result.matches) == len(embeddings)

    # Ensure `match` was called N times (depending on the input size as this is executed in batches under retryer)
    assert mock_matching_client.match.await_count == math.ceil(len(embeddings) / settings.GET_MATCHES_CHUNK_SIZE)


@pytest.mark.asyncio
async def test_get_knn_matches(
    populate_db: None,
    patch_matching_client: None,
    patch_bigquery_client: None,
    cell_operations_service_with_mocks: CellOperationsService,
    mock_valid_anndata: anndata.AnnData,
    cas_model: models.CASModel,
):
    """
    Test the `get_knn_matches` method.

    Validates that the method correctly retrieves nearest neighbors for all observations in the input AnnData object
    using the mocked `MatchingClient`.

    :param populate_db: Fixture to populate the test database.
    :param patch_matching_client: Fixture to patch the MatchingClient.
    :param patch_bigquery_client: Fixture to patch BigQuery client.
    :param cell_operations_service_with_mocks: Mocked CellOperationsService instance.
    :param mock_valid_anndata: Mocked AnnData object.
    :param cas_model: CASModel instance from the test database.
    """
    service = cell_operations_service_with_mocks
    query_ids, knn_response = await service.get_knn_matches(adata=mock_valid_anndata, model=cas_model)
    assert len(knn_response.matches) == len(mock_valid_anndata.obs)


@pytest.mark.asyncio
async def test_annotate_cell_type_summary_statistics_strategy_with_activity_logging_for_user_with_quota(
    populate_db: None,
    patch_matching_client: None,
    patch_bigquery_client: None,
    patch_starlette_context: None,
    db_session: Session,
    mock_valid_anndata: anndata.AnnData,
    mock_file_with_anndata: io.BytesIO,
    cell_operations_service_with_mocks: CellOperationsService,
    cas_model: models.CASModel,
    user_with_quota: models.User,
):
    """
    Tests the `annotate_cell_type_summary_statistics_strategy_with_activity_logging` method for a user with sufficient
    quota.

    Ensure that:
    - The response length matches the number of cells in the input AnnData object.
    - The `UserActivity` records are correctly created in the database, including both `STARTED` and `SUCCEEDED` events.

    :param populate_db: Use a fixture to populate the test database with mock data.
    :param patch_matching_client: Use a fixture to patch the `MatchingClient`.
    :param patch_bigquery_client: Use a fixture to patch the BigQuery client.
    :param patch_starlette_context: Use a fixture to patch the `starlette_context`.
    :param db_session: Database session fixture.
    :param mock_valid_anndata: Mocked AnnData object containing input data.
    :param mock_file_with_anndata: Mocked serialized AnnData file as input.
    :param cell_operations_service_with_mocks: Mocked `CellOperationsService` instance.
    :param cas_model: CASModel instance retrieved from the test database.
    :param user_with_quota: User instance with sufficient quota from the test database.
    """
    service = cell_operations_service_with_mocks

    cas_response = await service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        user=user_with_quota, file=mock_file_with_anndata, model_name=cas_model.model_name
    )
    # Check if the response length is the same as length of the input anndata file
    assert len(cas_response) == len(mock_valid_anndata)

    # Verify that UserActivity was created after the request
    user_activities = (
        db_session.query(models.UserActivity)
        .filter_by(
            user_id=user_with_quota.id,
            request_id=constants.TEST_REQUEST_ID_NEW,
            model_name=cas_model.model_name,
            method=constants.SUMMARY_STATS_METHOD,
        )
        .all()
    )
    assert len(user_activities) == 2
    first_user_activity = user_activities[0]
    second_user_activity = user_activities[1]
    assert first_user_activity.request_id == constants.TEST_REQUEST_ID_NEW
    assert second_user_activity.request_id == constants.TEST_REQUEST_ID_NEW
    assert first_user_activity.event == models.UserActivityEvent.STARTED
    assert second_user_activity.event == models.UserActivityEvent.SUCCEEDED


@pytest.mark.asyncio
async def test_annotate_cell_type_summary_statistics_strategy_with_activity_logging_for_user_without_quota(
    populate_db: None,
    patch_matching_client: None,
    patch_bigquery_client: None,
    mock_valid_anndata: anndata.AnnData,
    mock_file_with_anndata: io.BytesIO,
    cell_operations_service_with_mocks: CellOperationsService,
    db_session: Session,
    cas_model: models.CASModel,
    user_without_quota: models.User,
):
    """
    Tests the `annotate_cell_type_summary_statistics_strategy_with_activity_logging` method for a user without sufficient quota.

    Ensure that:
    - The method raises a `QuotaExceededException` when the user's quota is exceeded.
    - No `UserActivity` records are created in the database after the request.

    :param populate_db: Use a fixture to populate the test database with mock data.
    :param patch_matching_client: Use a fixture to patch the `MatchingClient`.
    :param patch_bigquery_client: Use a fixture to patch the BigQuery client.
    :param cell_operations_service_with_mocks: Mocked `CellOperationsService` instance.
    :param mock_valid_anndata: Mocked AnnData object containing input data.
    :param mock_file_with_anndata: Mocked serialized AnnData file as input.
    :param db_session: Database session fixture.
    :param cas_model: CASModel instance retrieved from the test database.
    :param user_without_quota: User instance without sufficient quota from the test database.
    """
    service = cell_operations_service_with_mocks
    # Use pytest.raises to check for QuotaExceededException
    with pytest.raises(exceptions.QuotaExceededException, match="User quota exceeded"):
        await service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
            user=user_without_quota, file=mock_file_with_anndata, model_name=cas_model.model_name
        )

    # Verify that UserActivity was not created after the request
    user_activities = (
        db_session.query(models.UserActivity)
        .filter_by(
            user_id=user_without_quota.id,
            request_id=constants.TEST_REQUEST_ID_NEW,
            model_name=cas_model.model_name,
            method=constants.SUMMARY_STATS_METHOD,
        )
        .all()
    )
    assert len(user_activities) == 0


@pytest.fixture
def patch_strategy_init_with_resource(
    cell_ontology_resource_mock: CellOntologyResource,
) -> t.Generator[Mock, None, None]:
    """
    Fixture to patch the `CellTypeOntologyAwareConsensusStrategy.__init__` method and inject a mocked resource.

    Allows overriding the initialization of the `CellTypeOntologyAwareConsensusStrategy` class to include the mocked
    ontology resource and other dependencies during testing.

    :param cell_ontology_resource_mock: Mocked ontology resource for the consensus strategy.
    """
    patch_ref = "casp.services.api.services.consensus_engine.CellTypeOntologyAwareConsensusStrategy.__init__"

    with patch(patch_ref, autospec=True) as mock_init:
        # Define a side effect that injects the mocked resource
        def init_side_effect(
            self,
            prune_threshold: float,
            weighting_prefactor: float,
            cell_ontology_resource=None,
            cell_operations_dm=None,
        ):
            self.prune_threshold = prune_threshold
            self.weighting_prefactor = weighting_prefactor
            self.cell_ontology_resource = cell_ontology_resource_mock
            self.cell_operations_dm = cell_operations_dm

        mock_init.side_effect = init_side_effect
        yield mock_init  # Return the patched constructor for assertions if needed


@pytest.mark.asyncio
async def test_annotate_cell_type_ontology_aware_strategy_with_activity_logging(
    populate_db: None,
    patch_matching_client: None,
    patch_bigquery_client: None,
    patch_starlette_context: None,
    patch_strategy_init_with_resource: Mock,
    mock_valid_anndata: anndata.AnnData,
    db_session: Session,
    mock_file_with_anndata: io.BytesIO,
    cell_operations_service_with_mocks: CellOperationsService,
    cell_ontology_resource_mock: CellOntologyResource,
    cas_model: models.CASModel,
    user_with_quota: models.User,
) -> None:
    """
    Tests the `annotate_cell_type_ontology_aware_strategy_with_activity_logging` method.

    Ensure that:
    - The response length matches the number of cells in the input AnnData object.
    - The `UserActivity` records are correctly created in the database, including both `STARTED` and `SUCCEEDED` events.
    - The consensus strategy is initialized correctly with the mocked ontology resource.

    :param populate_db: Use a fixture to populate the test database with mock data.
    :param patch_matching_client: Use a fixture to patch the `MatchingClient`.
    :param patch_bigquery_client: Use a fixture to patch the BigQuery client.
    :param patch_starlette_context: Use a fixture to patch the `starlette_context`.
    :param patch_strategy_init_with_resource: Use a fixture to patch the consensus strategy initialization.
    :param cell_operations_service_with_mocks: Mocked `CellOperationsService` instance.
    :param mock_valid_anndata: Mocked AnnData object containing input data.
    :param mock_file_with_anndata: Mocked serialized AnnData file as input.
    :param cell_ontology_resource_mock: Mocked ontology resource injected into the consensus strategy.
    :param cas_model: CASModel instance retrieved from the test database.
    :param user_with_quota: User instance with sufficient quota from the test database.
    :param db_session: Database session fixture.
    """
    service = cell_operations_service_with_mocks

    response = await service.annotate_cell_type_ontology_aware_strategy_with_activity_logging(
        user=user_with_quota,
        file=mock_file_with_anndata,
        model_name=cas_model.model_name,
        prune_threshold=0.1,
        weighting_prefactor=1.0,
    )

    # Validate the response has been generated correctly
    assert len(response) == len(mock_valid_anndata.obs)
    # Verify that UserActivity was created after the request
    user_activities = (
        db_session.query(models.UserActivity)
        .filter_by(
            user_id=user_with_quota.id,
            request_id=constants.TEST_REQUEST_ID_NEW,
            model_name=cas_model.model_name,
            method=constants.ONTOLOGY_AWARE_METHOD,
        )
        .all()
    )
    assert len(user_activities) == 2
    first_user_activity = user_activities[0]
    second_user_activity = user_activities[1]
    assert first_user_activity.request_id == constants.TEST_REQUEST_ID_NEW
    assert second_user_activity.request_id == constants.TEST_REQUEST_ID_NEW
    assert first_user_activity.event == models.UserActivityEvent.STARTED
    assert second_user_activity.event == models.UserActivityEvent.SUCCEEDED
    patch_strategy_init_with_resource.assert_called_once()
