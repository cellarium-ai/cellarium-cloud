"""
Tests things in the cell_operations_service that can reasonably be tested with a unit test.
"""

from unittest.mock import patch

import numpy as np
import pytest

from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.services import exceptions
from casp.services.api.services.cell_operations_service import CellOperationsService
from casp.services.constants import ContextKeys
from casp.services.db import models
from tests.unit import constants


@pytest.fixture(autouse=True)
def patch_matching_client(mock_matching_client):
    patching_method = "casp.services.api.clients.matching_client.MatchingClient.from_index"

    with patch(target=patching_method, return_value=mock_matching_client):
        yield


@pytest.fixture(autouse=True)
def mock_starlette_context():
    """
    Fixture to mock starlette_context.context globally.
    """
    from starlette_context.ctx import _request_scope_context_storage

    # Replace the internal context storage to avoid LookupError
    mock_data = {ContextKeys.sentry_trace_id: constants.TEST_REQUEST_ID_NEW}
    token = _request_scope_context_storage.set(mock_data)

    try:
        yield
    finally:
        _request_scope_context_storage.reset(token)


@pytest.fixture
def cas_model(db_session) -> models.CASModel:
    """
    Fixture to provide the CASModel instance from the database.

    :param db_session: Fixture providing the database session.
    :return: CASModel instance matching TEST_MODEL_NAME.
    """
    return db_session.query(models.CASModel).filter(models.CASModel.model_name == constants.TEST_MODEL_NAME).first()


@pytest.fixture
def user_with_quota(db_session) -> models.User:
    """
    Fixture to provide a user with sufficient quota from the database.

    :param db_session: Fixture providing the database session.
    :return: User instance with sufficient quota.
    """
    return db_session.query(models.User).filter(models.User.username == constants.USER_EMAIL_WITH_QUOTA).first()


@pytest.fixture
def user_without_quota(db_session) -> models.User:
    """
    Fixture to provide a user without sufficient quota from the database.

    :param db_session: Fixture providing the database session.
    :return: User instance without sufficient quota.
    """
    return db_session.query(models.User).filter(models.User.username == constants.USER_EMAIL_WITHOUT_QUOTA).first()


@pytest.fixture(
    params=[
        (np.random.rand(10, 3), 3, 4),  # 10 embeddings, chunk size 3
        (np.random.rand(15, 128), 5, 3),  # 15 embeddings, chunk size 5
        (np.random.rand(10, 64), 10, 1),  # Exactly one chunk
        (np.random.rand(10, 64), 20, 1),  # Chunk size larger than data
        (np.random.rand(0, 64), 5, 0),  # Empty embeddings
        (np.random.rand(11, 64), 4, 3),  # Irregular chunks
    ]
)
def embedding_test_data(request):
    """
    Fixture to provide test data for chunk splitting tests.

    :param request: Built-in pytest fixture to iterate through parameters.
    :return: Tuple of (embeddings, chunk_size, expected_num_chunks).
    """
    return request.param


def test__split_embeddings_into_chunks(embedding_test_data):
    """
    Test the __split_embeddings_into_chunks static method with various inputs.
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
    mock_matching_client, cell_operations_service_with_mocks, embedding_test_data
):
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
    mock_matching_client, cell_operations_service_with_mocks, db_session, cas_model
):
    embeddings = np.random.rand(10, 64)  # A chunk of 5 embeddings, each of 64 dimensions
    service = cell_operations_service_with_mocks
    result = await service.get_knn_matches_from_embeddings(embeddings=embeddings, model=cas_model)

    assert isinstance(result, MatchResult)
    assert len(result.matches) == len(embeddings)

    # Ensure `match` was called N times (depending on the input size as this is executed in batches under retryer
    assert mock_matching_client.match.await_count == 2


@pytest.mark.asyncio
async def test_get_knn_matches(
    mock_matching_client, cell_operations_service_with_mocks, mock_valid_anndata, db_session, cas_model
):
    service = cell_operations_service_with_mocks
    query_ids, knn_response = await service.get_knn_matches(adata=mock_valid_anndata, model=cas_model)
    assert len(knn_response.matches) == len(mock_valid_anndata.obs)


@pytest.mark.asyncio
async def test_annotate_cell_type_summary_statistics_strategy_with_activity_logging_for_user_with_quota(
    mock_matching_client,
    cell_operations_service_with_mocks,
    mock_valid_anndata,
    mock_file_with_anndata,
    mock_starlette_context,
    db_session,
    cas_model,
    user_with_quota,
):
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
    mock_matching_client,
    cell_operations_service_with_mocks,
    mock_valid_anndata,
    mock_file_with_anndata,
    db_session,
    cas_model,
    user_without_quota,
):
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
def mock_strategy_with_resource(cell_ontology_resource_mock):
    """
    Fixture to patch the CellTypeOntologyAwareConsensusStrategy.__init__ and inject the mocked resource.
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
    mock_matching_client,
    mock_strategy_with_resource,
    cell_operations_service_with_mocks,
    mock_valid_anndata,
    mock_file_with_anndata,
    cas_model,
    user_with_quota,
    cell_ontology_resource_mock,
    db_session,
):
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
    mock_strategy_with_resource.assert_called_once()
