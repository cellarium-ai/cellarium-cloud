from typing import NamedTuple, List
from casp.services.api import clients
from casp.services.db import models

from google.cloud import aiplatform_v1
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor


"""
Contains a common interface for querying the matching engine regarless of API (e.g. gRPC vs REST)
"""


class MatchResult(NamedTuple):

    class Neighbor(NamedTuple):
        """
        Neighbor represents a neighbor of a feature vector in the matching service.
        """
        id: str
        distance: float
        feature_vector: List[float] = None

    class NearestNeighbors(NamedTuple):
        """
        NearestNeighbors represents the result of a nearest neighbors request to the matching service. There should be 1
        per query that was sent.
        """
        neighbors: List["MatchResult.Neighbor"] = []

    """
    MatchResult represents the result of a match request to the matching service.
    """
    matches: List["MatchResult.NearestNeighbors"] = []

    def concat(self, other: "MatchResult") -> "MatchResult":
        """
        Concatenate the results of two match requests.

        :param other: The other match result to concatenate.

        :return: The concatenated match result as a new instance of MatchResult.
        """
        return MatchResult(matches=self.matches + other.matches)


class MatchingClient():
    def __init__(self, index: models.CASMatchingEngineIndex) -> None:
        self.index = index

    def match(self, queries: List[List[float]]) -> MatchResult:
        """
        Match queries against the specified index.

        :param index: The index to match against.
        :param queries: The queries to match.

        :return: The match result.
        """
        pass

    @staticmethod
    def from_index(index: models.CASMatchingEngineIndex) -> "MatchingClient":
        """
        Create a new MatchingClient from the specified index.

        :param index: The index to create the client from.

        :return: The new client.
        """
        if index.is_grpc:
            return MatchingClientGRPC(index=index)
        else:
            return MatchingClientREST(index=index)


class MatchingClientGRPC(MatchingClient):
    def __init__(self, index: models.CASMatchingEngineIndex) -> None:
        super().__init__(index)

    def match(self, queries: List[List[float]]) -> MatchResult:
        """
        Match queries against the specified index using gRPC.

        :param index: The index to match against.
        :param queries: The queries to match.

        :return: The match result.
        """
        index_endpoint_client = self.__get_match_index_endpoint_client()

        matches = index_endpoint_client.match(
            deployed_index_id=self.index.deployed_index_id,
            queries=queries,
            num_neighbors=self.index.num_neighbors,
        )

        return self.__adapt_result(matches)

    def __get_match_index_endpoint_client(
        self
    ) -> clients.CustomMatchingEngineIndexEndpointClient:
        return clients.CustomMatchingEngineIndexEndpointClient(index_endpoint_name=self.index.endpoint_id)

    def __adapt_result(self, result: List[List[MatchNeighbor]]) -> MatchResult:
        return MatchResult(
            matches=list(map(lambda r: MatchResult.NearestNeighbors(
                neighbors=list(map(lambda n: MatchResult.Neighbor(
                    id=n.id,
                    distance=n.distance,
                    feature_vector=n.feature_vector),
                    r))),
                result)))


class MatchingClientREST(MatchingClient):
    def __init__(self, index: models.CASMatchingEngineIndex) -> None:
        super().__init__(index)

    def match(self, queries: List[List[float]]) -> MatchResult:
        """
        Match queries against the specified index using REST.

        :param index: The index to match against.
        :param queries: The queries to match.

        :return: The match result.
        """

        vector_search_client = self.__get_match_index_endpoint_client()

        # Build query objects
        query_objects = list(map(lambda e: aiplatform_v1.FindNeighborsRequest.Query(
            datapoint=aiplatform_v1.IndexDatapoint(feature_vector=e), neighbor_count=self.index.num_neighbors), queries))

        # Prepare the request to be sent
        request = aiplatform_v1.FindNeighborsRequest(
            index_endpoint=self.index.endpoint_id,
            deployed_index_id=self.index.deployed_index_id,
            queries=query_objects,
            return_full_datapoint=True,
        )

        # Send the request
        matches = vector_search_client.find_neighbors(request)

        return self.__adapt_result(matches)

    # @classmethod
    def __get_match_index_endpoint_client(
        self
    ) -> aiplatform_v1.MatchServiceClient:
        client_options = {
            "api_endpoint": self.index.api_endpoint
        }
        return aiplatform_v1.MatchServiceClient(
            client_options=client_options,
        )

    def __adapt_result(self, result: aiplatform_v1.FindNeighborsResponse) -> MatchResult:
        return MatchResult(
            matches=list(map(lambda q: MatchResult.NearestNeighbors(
                neighbors=list(map(lambda n: MatchResult.Neighbor(
                    id=n.datapoint.datapoint_id,
                    distance=n.distance,
                    feature_vector=n.datapoint.feature_vector),
                    q.neighbors))),
                result.nearest_neighbors)))
