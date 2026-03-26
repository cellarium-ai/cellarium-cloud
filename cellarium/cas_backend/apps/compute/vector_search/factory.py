from cellarium.cas_backend.apps.compute.vector_search.exceptions import VectorSearchConfigurationError
from cellarium.cas_backend.apps.compute.vector_search.protocol import VectorSearchProtocol
from cellarium.cas_backend.apps.compute.vector_search.tiledb import TileDBVectorSearch
from cellarium.cas_backend.apps.compute.vector_search.vertex_ai import (
    VertexVectorSearchClientGRPC,
    VertexVectorSearchClientREST,
)
from cellarium.cas_backend.core.db import models


def from_index(index: models.CASVectorIndex | models.CASMatchingEngineIndex) -> VectorSearchProtocol:
    if isinstance(index, models.CASVectorIndex):
        return TileDBVectorSearch(index=index)

    if isinstance(index, models.CASMatchingEngineIndex):
        if index.is_grpc:
            return VertexVectorSearchClientGRPC(index=index)
        return VertexVectorSearchClientREST(index=index)

    raise VectorSearchConfigurationError("Unsupported vector search index configuration.")


def from_model(model: models.CASModel) -> VectorSearchProtocol:
    if getattr(model, "cas_vector_index", None) is not None:
        return from_index(model.cas_vector_index)

    if getattr(model, "cas_matching_engine", None) is not None:
        return from_index(model.cas_matching_engine)

    raise VectorSearchConfigurationError(f"Model '{model.model_name}' does not have any vector search configuration.")
