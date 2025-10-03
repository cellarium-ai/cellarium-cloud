import asyncio

from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.services.annotation_engines.ann_consensus_engine import ConsensusStrategyInterface


class ConsensusEngine:
    """
    Consensus engine responsible for summarizing query neighbor context using a specified strategy.
    """

    def __init__(self, strategy: ConsensusStrategyInterface):
        self.strategy = strategy

    def summarize(self, knn_query: MatchResult) -> schemas.QueryAnnotationAbstractType:
        """
        Summarize query neighbor context using the specified strategy.

        :param knn_query: The result of the kNN query.

        :return: A list of summarized query neighbor context annotations.
        """
        return self.strategy.summarize(knn_query=knn_query)

    async def summarize_async(self, knn_query: MatchResult) -> schemas.QueryAnnotationAbstractType:
        """
        Call :meth:`summarize` in a separate thread to avoid blocking the event loop.

        :param knn_query: The result of the kNN query.

        :return: A list of summarized query neighbor context annotations.
        """
        return await asyncio.to_thread(self.summarize, knn_query)
