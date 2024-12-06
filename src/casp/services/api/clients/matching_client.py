import typing as t

from google.cloud import aiplatform_v1
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor
from pydantic import BaseModel, Field

from casp.services.api import clients
from casp.services.db import models


class MatchResult(BaseModel):
    """
    MatchResult represents the result of a match request to the matching service.
    """

    class Neighbor(BaseModel):
        """
        Neighbor represents a neighbor of a feature vector (e.g. a cell) in the matching service.
        """

        cas_cell_index: str = Field(default=None)
        distance: float = Field(default=None)

    class NearestNeighbors(BaseModel):
        """
        NearestNeighbors represents the result of a nearest neighbors request to the matching service. There should be 1
        per query that was sent.
        """

        neighbors: t.List["MatchResult.Neighbor"] = Field(default=[])

    matches: t.List["MatchResult.NearestNeighbors"] = Field(default=[])

    def concat(self, other: "MatchResult") -> "MatchResult":
        """
        Concatenate the results of two match requests.

        :param other: The other match result to concatenate.

        :return: The concatenated match result as a new instance of MatchResult.
        """
        return MatchResult(matches=self.matches + other.matches)

    def get_unique_ids(self) -> t.List[int]:
        """
        Get unique cas cell ids across all neighbors in querying cells in :class:`MatchResult`

        :return: A list with unique indexes of neighbors
        """
        return list(
            set(
                int(neighbor.cas_cell_index)
                for query_neighbors in self.matches
                for neighbor in query_neighbors.neighbors
            )
        )


MatchResult.NearestNeighbors.update_forward_refs()
MatchResult.Neighbor.update_forward_refs()


class MatchingClient:
    """
    Common interface for querying the matching service regardless of API (e.g. gRPC vs REST)
    """

    def __init__(self, index: models.CASMatchingEngineIndex):
        self.index = index

    async def match(self, queries: t.List[t.List[float]]) -> MatchResult:
        """
        Match queries against the specified index.

        :param queries: The queries to match.

        :return: The match result.
        """
        raise NotImplementedError

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
    """
    MatchingClientGRPC is a client for matching queries against a Vertex AI vector index when the index
    has been deployed in private mode, which happens to use a gRPC client.
    """

    def __init__(self, index: models.CASMatchingEngineIndex):
        super().__init__(index=index)
        self.index_endpoint_client = clients.CustomMatchingEngineIndexEndpointClient(
            index_endpoint_name=self.index.endpoint_id
        )

    def __adapt_result(self, result: t.List[t.List[MatchNeighbor]]) -> MatchResult:
        matches = []
        for query_result in result:
            neighbors = []
            for neighbor in query_result:
                # TODO: We have to make sure that index returns a distance, not a similarity score. All indexes
                # that we use at the moment (DOT_PRODUCT_DISTANCE with L2_NORM_TYPE) return cosine similarity instead
                # of cosine distance. We have to develop a more sophisticated handler for index depending on its type.
                # Also we should get rid of this boilerplate code and comments..
                distance = round(1 - neighbor.distance, 9)
                neighbors.append(
                    MatchResult.Neighbor(
                        cas_cell_index=neighbor.id,
                        distance=distance,
                    )
                )
            matches.append(MatchResult.NearestNeighbors(neighbors=neighbors))
        return MatchResult(matches=matches)

    async def match(self, queries: t.List[t.List[float]]) -> MatchResult:
        """
        Match queries against the specified index using gRPC.

        :param queries: The queries to match.

        :return: The match result.
        """
        matches = self.index_endpoint_client.match(
            deployed_index_id=self.index.deployed_index_id, queries=queries, num_neighbors=self.index.num_neighbors
        )
        return self.__adapt_result(matches)


class MatchingClientREST(MatchingClient):
    """
    MatchingClientREST is a client for matching queries against a Vertex AI vector index when the index
    has been deployed in public mode, which happens to use a REST client.
    """

    def __init__(self, index: models.CASMatchingEngineIndex):
        super().__init__(index=index)
        client_options = {"api_endpoint": self.index.api_endpoint}
        self.vector_search_client = aiplatform_v1.MatchServiceAsyncClient(client_options=client_options)

    def __adapt_result(self, result: aiplatform_v1.FindNeighborsResponse) -> MatchResult:
        matches = []
        for query_reqult in result.nearest_neighbors:
            neighbors = []
            for neighbor in query_reqult.neighbors:
                # TODO: We have to make sure that index returns a distance, not a similarity score. All indexes
                # that we use at the moment (DOT_PRODUCT_DISTANCE with L2_NORM_TYPE return cosine similarity instead of
                # cosine distance. We have to develop a more sophisticated handler for index depending on its type.
                # Also we should get rid of this boilerplate code and comments.
                distance = round(1 - neighbor.distance, 9)
                neighbors.append(
                    MatchResult.Neighbor(
                        cas_cell_index=neighbor.datapoint.datapoint_id,
                        distance=distance,
                    )
                )
            matches.append(MatchResult.NearestNeighbors(neighbors=neighbors))
        return MatchResult(matches=matches)

    async def match(self, queries: t.List[t.List[float]]) -> MatchResult:
        """
        Match queries against the specified index using REST.

        :param queries: The queries to match.

        :return: The match result.
        """
        # Build query objects
        query_objects = [
            aiplatform_v1.FindNeighborsRequest.Query(
                datapoint=aiplatform_v1.IndexDatapoint(feature_vector=e), neighbor_count=self.index.num_neighbors
            )
            for e in queries
        ]

        # Prepare the request to be sent
        request = aiplatform_v1.FindNeighborsRequest(
            index_endpoint=self.index.endpoint_id,
            deployed_index_id=self.index.deployed_index_id,
            queries=query_objects,
            return_full_datapoint=False,
        )

        # Send the request
        matches = await self.vector_search_client.find_neighbors(request)
        return self.__adapt_result(matches)
