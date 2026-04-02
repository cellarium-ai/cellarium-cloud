from cellarium.cas_backend.apps.compute.vector_search.exceptions import VectorSearchConfigurationError
from cellarium.cas_backend.apps.compute.vector_search.tiledb import TileDBVectorSearch
from cellarium.cas_backend.core.db import models


def from_model(model: models.CASModel) -> TileDBVectorSearch:
    if getattr(model, "cas_vector_index", None) is not None:
        return TileDBVectorSearch(index=model.cas_vector_index)

    raise VectorSearchConfigurationError(f"Model '{model.model_name}' does not have any vector search configuration.")
