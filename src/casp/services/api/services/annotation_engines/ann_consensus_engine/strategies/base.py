import typing as t

from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult


class ConsensusStrategyInterface:
    def summarize(self, knn_query: MatchResult) -> t.List[schemas.QueryCellAnnotationAbstract]:
        """
        Define the method signature for the consensus strategy, which is responsible for summarizing query neighbor

        :param knn_query: The result of the kNN query.

        :return: A list of summarized query neighbor context annotations.
        """
        raise NotImplementedError
