from google.cloud import aiplatform_v1
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor
import numpy as np

from cellarium.cas_backend.apps.compute.vector_search.protocol import VectorSearchProtocol
from cellarium.cas_backend.apps.compute.vector_search.schemas import MatchResult
from cellarium.cas_backend.apps.compute.vector_search.vertex_ai import _matching_engine
from cellarium.cas_backend.core.db import models


class VertexVectorSearchClientGRPC(VectorSearchProtocol):
    """
    Client for querying a Vertex AI vector index deployed in private mode via gRPC.
    """

    def __init__(self, index: models.CASMatchingEngineIndex):
        self.index = index
        self.index_endpoint_client = _matching_engine.CustomMatchingEngineIndexEndpointClient(
            index_endpoint_name=self.index.endpoint_id
        )

    def __adapt_result(self, result: list[list[MatchNeighbor]]) -> MatchResult:
        matches = []
        for query_result in result:
            neighbors = []
            for neighbor in query_result:
                distance = round(1 - neighbor.distance, 9)
                neighbors.append(
                    MatchResult.Neighbor(
                        cas_cell_index=neighbor.id,
                        distance=distance,
                    )
                )
            matches.append(MatchResult.NearestNeighbors(neighbors=neighbors))
        return MatchResult(matches=matches)

    async def match(self, embeddings: np.ndarray) -> MatchResult:
        matches = self.index_endpoint_client.match(
            deployed_index_id=self.index.deployed_index_id,
            queries=embeddings.tolist(),
            num_neighbors=self.index.num_neighbors,
        )
        return self.__adapt_result(matches)


class VertexVectorSearchClientREST(VectorSearchProtocol):
    """
    Client for querying a Vertex AI vector index deployed in public mode via REST.
    """

    def __init__(self, index: models.CASMatchingEngineIndex):
        self.index = index
        client_options = {"api_endpoint": self.index.api_endpoint}
        self.vector_search_client = aiplatform_v1.MatchServiceAsyncClient(client_options=client_options)

    def __adapt_result(self, result: aiplatform_v1.FindNeighborsResponse) -> MatchResult:
        matches = []
        for query_result in result.nearest_neighbors:
            neighbors = []
            for neighbor in query_result.neighbors:
                distance = round(1 - neighbor.distance, 9)
                neighbors.append(
                    MatchResult.Neighbor(
                        cas_cell_index=neighbor.datapoint.datapoint_id,
                        distance=distance,
                    )
                )
            matches.append(MatchResult.NearestNeighbors(neighbors=neighbors))
        return MatchResult(matches=matches)

    async def match(self, embeddings: np.ndarray) -> MatchResult:
        query_objects = [
            aiplatform_v1.FindNeighborsRequest.Query(
                datapoint=aiplatform_v1.IndexDatapoint(feature_vector=e), neighbor_count=self.index.num_neighbors
            )
            for e in embeddings.tolist()
        ]

        request = aiplatform_v1.FindNeighborsRequest(
            index_endpoint=self.index.endpoint_id,
            deployed_index_id=self.index.deployed_index_id,
            queries=query_objects,
            return_full_datapoint=False,
        )

        matches = await self.vector_search_client.find_neighbors(request)
        return self.__adapt_result(matches)
