from typing import Protocol

from cellarium.cas_backend.apps.compute.vector_search.schemas import MatchResult


class VectorSearchProtocol(Protocol):
    """
    Common protocol for querying the vector search backend.
    """

    async def match(self, queries: list[list[float]]) -> MatchResult:
        """
        Match queries against the configured vector search backend.

        :param queries: The queries to match.

        :return: The match result.
        """
