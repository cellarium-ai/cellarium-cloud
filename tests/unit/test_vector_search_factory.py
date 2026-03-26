from unittest.mock import Mock, patch

import pytest

from cellarium.cas_backend.apps.compute.vector_search.exceptions import VectorSearchConfigurationError
from cellarium.cas_backend.apps.compute.vector_search.factory import from_model
from cellarium.cas_backend.core.db import models
from tests.unit.fixtures import constants


def _build_model() -> models.CASModel:
    return models.CASModel(
        id=1,
        model_name=constants.TEST_MODEL_NAME,
        model_file_path=constants.TEST_MODEL_FILE_PATH,
        embedding_dimension=constants.TEST_EMBEDDING_DIMENSION,
    )


def _build_matching_engine() -> models.CASMatchingEngineIndex:
    return models.CASMatchingEngineIndex(
        id=1,
        index_name=constants.TEST_INDEX_NAME,
        model_id=1,
        num_neighbors=constants.TEST_INDEX_NUM_NEIGHBORS,
        embedding_dimension=constants.TEST_EMBEDDING_DIMENSION,
        endpoint_id=constants.TEST_INDEX_ENDPOINT_ID,
        is_grpc=True,
    )


def _build_vector_index(**overrides) -> models.CASVectorIndex:
    payload = {
        "id": 1,
        "index_name": constants.TEST_VECTOR_INDEX_NAME,
        "model_id": 1,
        "num_neighbors": constants.TEST_INDEX_NUM_NEIGHBORS,
        "embedding_dimension": constants.TEST_EMBEDDING_DIMENSION,
        "index_uri": constants.TEST_VECTOR_INDEX_URI,
        "index_type": constants.TEST_VECTOR_INDEX_TYPE,
        "distance_metric": constants.TEST_VECTOR_DISTANCE_METRIC,
        "nprobe": constants.TEST_VECTOR_INDEX_NPROBE,
        "l_search": constants.TEST_VECTOR_INDEX_L_SEARCH,
    }
    payload.update(overrides)
    return models.CASVectorIndex(**payload)


def test_from_model_prefers_tiledb_index() -> None:
    model = _build_model()
    model.cas_vector_index = _build_vector_index()
    model.cas_matching_engine = _build_matching_engine()

    expected_client = Mock()
    with (
        patch(
            "cellarium.cas_backend.apps.compute.vector_search.factory.TileDBVectorSearch",
            return_value=expected_client,
        ) as tiledb_client,
        patch("cellarium.cas_backend.apps.compute.vector_search.factory.VertexVectorSearchClientGRPC") as grpc_client,
    ):
        assert from_model(model) is expected_client

    tiledb_client.assert_called_once_with(index=model.cas_vector_index)
    grpc_client.assert_not_called()


def test_from_model_uses_vertex_when_no_vector_index() -> None:
    model = _build_model()
    model.cas_matching_engine = _build_matching_engine()

    expected_client = Mock()
    with patch(
        "cellarium.cas_backend.apps.compute.vector_search.factory.VertexVectorSearchClientGRPC",
        return_value=expected_client,
    ) as grpc_client:
        assert from_model(model) is expected_client

    grpc_client.assert_called_once_with(index=model.cas_matching_engine)


def test_from_model_raises_when_no_index_exists() -> None:
    model = _build_model()

    with pytest.raises(VectorSearchConfigurationError, match="does not have any vector search"):
        from_model(model)


def test_from_model_invalid_vector_index_does_not_fallback() -> None:
    model = _build_model()
    model.cas_vector_index = _build_vector_index(distance_metric="invalid")
    model.cas_matching_engine = _build_matching_engine()

    with pytest.raises(VectorSearchConfigurationError, match="distance_metric"):
        from_model(model)
