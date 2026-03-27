from typing import Protocol

import numpy as np

from cellarium.cas_backend.apps.compute.vector_search.schemas import MatchResult


class VectorSearchProtocol(Protocol):
    """
    Common protocol for querying the vector search backend.
    """

    async def match(self, embeddings: np.ndarray) -> MatchResult:
        """
        Match embeddings against the configured vector search backend.

        :param embeddings: 2D float32 array of shape (n_cells, embedding_dim).

        :return: The match result.
        """
