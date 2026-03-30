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
        cell_metadata_uri=constants.TEST_MODEL_CELL_METADATA_URI,
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


def test_from_model_creates_tiledb_client() -> None:
    model = _build_model()
    model.cas_vector_index = _build_vector_index()

    expected_client = Mock()
    with patch(
        "cellarium.cas_backend.apps.compute.vector_search.factory.TileDBVectorSearch",
        return_value=expected_client,
    ) as tiledb_client:
        assert from_model(model) is expected_client

    tiledb_client.assert_called_once_with(index=model.cas_vector_index)


def test_from_model_raises_when_no_index_exists() -> None:
    model = _build_model()

    with pytest.raises(VectorSearchConfigurationError, match="does not have any vector search"):
        from_model(model)


def test_from_model_raises_for_invalid_vector_index() -> None:
    model = _build_model()
    model.cas_vector_index = _build_vector_index(distance_metric="invalid")

    with pytest.raises(VectorSearchConfigurationError, match="distance_metric"):
        from_model(model)
