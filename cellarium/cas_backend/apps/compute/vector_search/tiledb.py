from collections.abc import Mapping
from dataclasses import dataclass
import typing as t

import numpy as np
from tiledb import vector_search

from cellarium.cas_backend.apps.compute.vector_search.exceptions import VectorSearchConfigurationError
from cellarium.cas_backend.apps.compute.vector_search.schemas import MatchResult
from cellarium.cas_backend.core.db import models

INDEX_TYPE_FLAT = "FLAT"
INDEX_TYPE_IVF_FLAT = "IVF_FLAT"
INDEX_TYPE_VAMANA = "VAMANA"
SUPPORTED_INDEX_TYPES = {INDEX_TYPE_FLAT, INDEX_TYPE_IVF_FLAT, INDEX_TYPE_VAMANA}
SUPPORTED_DISTANCE_METRICS = {"cosine", "l2", "dot_product"}

_INDEX_CACHE: dict[str, t.Any] = {}


@dataclass(frozen=True)
class TileDBIndexConfig:
    index_uri: str
    index_type: str
    distance_metric: str
    nprobe: int | None
    l_search: int | None


def _get_index_class(index_type: str):
    match index_type:
        case "FLAT":
            index_class = vector_search.flat_index.FlatIndex
        case "IVF_FLAT":
            index_class = vector_search.ivf_flat_index.IVFFlatIndex
        case "VAMANA":
            index_class = vector_search.vamana_index.VamanaIndex
        case _:
            raise ValueError(f"Unsupported TileDB index_type '{index_type}'.")

    return index_class


def _validate_index_type(index_type: t.Any) -> str:
    if not isinstance(index_type, str):
        raise ValueError("TileDB vector index_type must be a string.")

    normalized = index_type.upper()
    if normalized not in SUPPORTED_INDEX_TYPES:
        raise ValueError(f"Unsupported TileDB index_type '{index_type}'.")
    return normalized


def validate_tiledb_index(index: models.CASVectorIndex) -> TileDBIndexConfig:
    if not isinstance(index.index_uri, str) or not index.index_uri.strip():
        raise ValueError("TileDB index_uri must be a non-empty string.")

    distance_metric = index.distance_metric
    if not isinstance(distance_metric, str) or distance_metric not in SUPPORTED_DISTANCE_METRICS:
        raise ValueError("TileDB distance_metric must be one of: cosine, l2, dot_product.")

    index_type = _validate_index_type(index.index_type)

    if index_type == INDEX_TYPE_IVF_FLAT:
        nprobe = index.nprobe
        if not isinstance(nprobe, int) or nprobe <= 0:
            raise ValueError("TileDB IVF_FLAT nprobe must be a positive integer.")
    if index_type == INDEX_TYPE_VAMANA:
        l_search = index.l_search
        if not isinstance(l_search, int) or l_search <= 0:
            raise ValueError("TileDB VAMANA l_search must be a positive integer.")

    return TileDBIndexConfig(
        index_uri=index.index_uri,
        index_type=index_type,
        distance_metric=distance_metric,
        nprobe=index.nprobe,
        l_search=index.l_search,
    )


def _get_cached_index(index_uri: str, index_type: str):
    if index_uri in _INDEX_CACHE:
        return _INDEX_CACHE[index_uri]

    index_class = _get_index_class(index_type)
    index_obj = index_class(uri=index_uri)

    _INDEX_CACHE[index_uri] = index_obj
    return index_obj


def _coerce_result_arrays(result: t.Any) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(result, Mapping):
        distances = result.get("distances")
        ids = result.get("ids")
    elif hasattr(result, "distances") and hasattr(result, "ids"):
        distances = result.distances
        ids = result.ids
    elif isinstance(result, tuple | list) and len(result) == 2:
        first, second = result
        first_array = np.asarray(first)
        second_array = np.asarray(second)

        if np.issubdtype(first_array.dtype, np.integer) and not np.issubdtype(second_array.dtype, np.integer):
            ids, distances = first_array, second_array
        else:
            distances, ids = first_array, second_array
    else:
        raise ValueError("Unexpected TileDB query result format.")

    distances_array = np.asarray(distances)
    ids_array = np.asarray(ids)
    if distances_array.ndim != 2 or ids_array.ndim != 2:
        raise ValueError("TileDB query results must be two-dimensional arrays.")
    if distances_array.shape != ids_array.shape:
        raise ValueError("TileDB query result ids and distances must have the same shape.")
    return distances_array, ids_array


class TileDBVectorSearch:
    """
    Client for querying TileDB vector indexes from the compute service.
    """

    def __init__(self, index: models.CASVectorIndex):
        self.index = index
        try:
            self.index_config = validate_tiledb_index(index)
        except (RuntimeError, ValueError) as exc:
            raise VectorSearchConfigurationError(str(exc)) from exc

        self.index_obj = _get_cached_index(
            index_uri=self.index_config.index_uri,
            index_type=self.index_config.index_type,
        )

        dimensions = self.index_obj.get_dimensions()
        if dimensions != self.index.embedding_dimension:
            raise VectorSearchConfigurationError(
                f"TileDB index dimension {dimensions} does not match configured embedding dimension "
                f"{self.index.embedding_dimension}."
            )

    def _query_sync(self, queries: np.ndarray):
        query_kwargs = {"k": self.index.num_neighbors}
        index_type = self.index_config.index_type
        if index_type == INDEX_TYPE_IVF_FLAT:
            query_kwargs["nprobe"] = self.index_config.nprobe
        elif index_type == INDEX_TYPE_VAMANA:
            query_kwargs["l_search"] = self.index_config.l_search

        return self.index_obj.query_internal(queries=queries, **query_kwargs)

    def _adapt_result(self, result: t.Any) -> MatchResult:
        try:
            distances, ids = _coerce_result_arrays(result)
        except ValueError as exc:
            raise VectorSearchConfigurationError(str(exc)) from exc

        matches = []
        for query_ids, query_distances in zip(ids, distances, strict=False):
            neighbors = [
                MatchResult.Neighbor(cas_cell_index=int(cas_cell_index), distance=round(float(distance), 9))
                for cas_cell_index, distance in zip(query_ids, query_distances, strict=False)
            ]
            matches.append(MatchResult.NearestNeighbors(neighbors=neighbors))
        return MatchResult(matches=matches)

    def match(self, embeddings: np.ndarray) -> MatchResult:
        if embeddings.ndim != 2:
            raise VectorSearchConfigurationError("TileDB embeddings must be a 2D array.")
        if embeddings.shape[1] != self.index.embedding_dimension:
            raise VectorSearchConfigurationError(
                f"Embedding dimension {embeddings.shape[1]} does not match configured embedding dimension "
                f"{self.index.embedding_dimension}."
            )

        result = self._query_sync(embeddings)
        return self._adapt_result(result)
