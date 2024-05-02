import asyncio
import typing as t

from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.services.consensus_engine.strategies import ConsensusStrategyProtocol


class ConsensusEngine:
    """
    Consensus engine responsible for summarizing query neighbor context using a specified strategy.
    """

    def __init__(self, strategy: ConsensusStrategyProtocol):
        self.strategy = strategy

    def summarize(self, query_ids: t.List[str], knn_query: MatchResult) -> schemas.QueryAnnotationAbstractType:
        """
        Summarize query neighbor context using the specified strategy.

        :param query_ids: List of query cell IDs.
        :param knn_query: The result of the kNN query.

        :return: A list of summarized query neighbor context annotations.
        """
        return self.strategy.summarize(query_ids, knn_query)

    async def summarize_async(
        self, query_ids: t.List[str], knn_query: MatchResult
    ) -> schemas.QueryAnnotationAbstractType:
        """
        Call :meth:`summarize` in a separate thread to avoid blocking the event loop.

        :param query_ids: List of query cell IDs.
        :param knn_query: The result of the kNN query.

        :return: A list of summarized query neighbor context annotations.
        """
        return await asyncio.to_thread(self.summarize, query_ids, knn_query)
