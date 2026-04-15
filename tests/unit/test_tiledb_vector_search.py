import numpy as np
import pytest

from cellarium.cas_backend.apps.compute.vector_search.exceptions import VectorSearchConfigurationError
from cellarium.cas_backend.apps.compute.vector_search.tiledb import (
    _INDEX_CACHE,
    TileDBVectorSearch,
)
from cellarium.cas_backend.core.db import models
from tests.unit.fixtures import constants


class FakeIVFFlatIndex:
    def __init__(self, uri, **kwargs):
        self.uri = uri
        self.init_kwargs = kwargs
        self.query_calls = []

    def get_dimensions(self):
        return constants.TEST_EMBEDDING_DIMENSION

    def query_internal(self, queries, **kwargs):
        self.query_calls.append((queries, kwargs))
        query_count = len(queries)
        return (
            np.array([[0.1, 0.2, 0.3]] * query_count, dtype=np.float32),
            np.array([[1, 2, 3]] * query_count, dtype=np.int64),
        )


class FakeDimensionMismatchIndex(FakeIVFFlatIndex):
    def get_dimensions(self):
        return constants.TEST_EMBEDDING_DIMENSION + 1


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
        "memory_budget": None,
    }
    payload.update(overrides)
    return models.CASVectorIndex(**payload)


@pytest.fixture(autouse=True)
def clear_tiledb_cache():
    _INDEX_CACHE.clear()
    yield
    _INDEX_CACHE.clear()


def test_init_raises_for_invalid_distance_metric() -> None:
    with pytest.raises(VectorSearchConfigurationError, match="distance_metric"):
        TileDBVectorSearch(index=_build_vector_index(distance_metric="bad-metric"))


def test_init_raises_for_dimension_mismatch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cellarium.cas_backend.apps.compute.vector_search.tiledb.vector_search.ivf_flat_index.IVFFlatIndex",
        FakeDimensionMismatchIndex,
    )

    with pytest.raises(
        VectorSearchConfigurationError,
        match="does not match configured embedding dimension",
    ):
        TileDBVectorSearch(index=_build_vector_index())


def test_match_adapts_multi_query_results(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cellarium.cas_backend.apps.compute.vector_search.tiledb.vector_search.ivf_flat_index.IVFFlatIndex",
        FakeIVFFlatIndex,
    )

    client = TileDBVectorSearch(index=_build_vector_index())
    embeddings = np.array(
        [[0.0] * constants.TEST_EMBEDDING_DIMENSION, [1.0] * constants.TEST_EMBEDDING_DIMENSION], dtype=np.float32
    )

    result = client.match(embeddings=embeddings)

    assert len(result.matches) == 2
    assert [neighbor.cas_cell_index for neighbor in result.matches[0].neighbors] == [1, 2, 3]
    assert [neighbor.distance for neighbor in result.matches[0].neighbors] == pytest.approx([0.1, 0.2, 0.3])
    assert client.index_obj.query_calls[0][1]["k"] == constants.TEST_INDEX_NUM_NEIGHBORS
    assert client.index_obj.query_calls[0][1]["nprobe"] == constants.TEST_VECTOR_INDEX_NPROBE


def test_match_raises_for_query_dimension_mismatch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cellarium.cas_backend.apps.compute.vector_search.tiledb.vector_search.ivf_flat_index.IVFFlatIndex",
        FakeIVFFlatIndex,
    )
    client = TileDBVectorSearch(index=_build_vector_index())

    with pytest.raises(VectorSearchConfigurationError, match="Embedding dimension"):
        client.match(embeddings=np.array([[1.0, 2.0, 3.0]], dtype=np.float32))


def test_memory_budget_none_does_not_pass_kwarg(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cellarium.cas_backend.apps.compute.vector_search.tiledb.vector_search.ivf_flat_index.IVFFlatIndex",
        FakeIVFFlatIndex,
    )

    client = TileDBVectorSearch(index=_build_vector_index(memory_budget=None))

    assert "memory_budget" not in client.index_obj.init_kwargs


def test_memory_budget_provided_passes_kwarg(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "cellarium.cas_backend.apps.compute.vector_search.tiledb.vector_search.ivf_flat_index.IVFFlatIndex",
        FakeIVFFlatIndex,
    )

    client = TileDBVectorSearch(index=_build_vector_index(memory_budget=4_000_000))

    assert client.index_obj.init_kwargs["memory_budget"] == 4_000_000
