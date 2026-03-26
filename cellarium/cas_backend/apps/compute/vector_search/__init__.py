from cellarium.cas_backend.apps.compute.vector_search.factory import from_index, from_model
from cellarium.cas_backend.apps.compute.vector_search.protocol import VectorSearchProtocol
from cellarium.cas_backend.apps.compute.vector_search.schemas import MatchResult
from cellarium.cas_backend.apps.compute.vector_search.tiledb import TileDBVectorSearch
from cellarium.cas_backend.apps.compute.vector_search.vertex_ai import (
    VertexVectorSearchClientGRPC,
    VertexVectorSearchClientREST,
)

__all__ = [
    "VectorSearchProtocol",
    "TileDBVectorSearch",
    "VertexVectorSearchClientGRPC",
    "VertexVectorSearchClientREST",
    "MatchResult",
    "from_index",
    "from_model",
]
