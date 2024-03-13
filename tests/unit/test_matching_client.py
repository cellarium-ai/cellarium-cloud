from typing import List
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor
from google.cloud.aiplatform_v1 import MatchServiceClient, FindNeighborsResponse, FindNeighborsRequest, IndexDatapoint
from parameterized import parameterized
from google.protobuf.json_format import Parse
from mockito import when, mock, unstub


from casp.services.db import models
from casp.services.api.clients import CustomMatchingEngineIndexEndpointClient
from casp.services.api.clients.matching import matching_client
from casp.services.api.clients.matching.matching_client import MatchingClient, MatchingClientGRPC, MatchingClientREST, MatchResult

import tests.unit.test_utils as test_utils


GRPC_INDEX = models.CASMatchingEngineIndex(
    is_grpc=True, endpoint_id="endpoint_id_grpc", deployed_index_id="deployed_index_id_grpc", num_neighbors=3)
REST_INDEX = models.CASMatchingEngineIndex(
    is_grpc=False, endpoint_id="endpoint_id_rest", deployed_index_id="deployed_index_id_rest", api_endpoint="api_endpoint_rest", num_neighbors=3)


class TestMatchingClient:
    """
    Test the MatchingClient class.

    """

    def teardown_method(self) -> None:
        unstub()

    def test_client_selector(self) -> None:
        """
        Test the creation of a matching clients from gRPC and REST indexes.

        """
        isinstance(MatchingClient.from_index(GRPC_INDEX), MatchingClientGRPC)
        isinstance(MatchingClient.from_index(REST_INDEX), MatchingClientREST)

    def test_concat_matches(self) -> None:
        """
        Test concatenating match results.

        """
        match1 = MatchResult(matches=[MatchResult.NearestNeighbors(
            neighbors=[MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2])])])
        match2 = MatchResult(matches=[MatchResult.NearestNeighbors(
            neighbors=[MatchResult.Neighbor(id="1", distance=1.0, feature_vector=[10, 11, 12])])])

        assert match1.concat(match2) == MatchResult(matches=[
            MatchResult.NearestNeighbors(neighbors=[
                MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2]),
            ]),
            MatchResult.NearestNeighbors(neighbors=[
                MatchResult.Neighbor(id="1", distance=1.0, feature_vector=[10, 11, 12]),
            ]),
        ])

    @parameterized.expand([
        ([], [], FindNeighborsResponse(), MatchResult()),
        (
            [[1, 2, 3]],
            [FindNeighborsRequest.Query(neighbor_count=3, datapoint=IndexDatapoint(feature_vector=[1, 2, 3]))],
            # Note: reading from JSON so that we have a few examples of responses from the REST service.
            Parse(test_utils.read_resource("tests/unit/test_query_responses/rest_response_0.json"), FindNeighborsResponse()._pb),
            MatchResult(matches=[
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2]),
                    matching_client.MatchResult.Neighbor(id="1", distance=10.0, feature_vector=[3, 4, 5]),
                    matching_client.MatchResult.Neighbor(id="2", distance=20.0, feature_vector=[6, 7, 8]),
                ])
            ])
        ),
        (
            [[1, 2, 3], [4, 5, 6]],
            [
                FindNeighborsRequest.Query(neighbor_count=3, datapoint=IndexDatapoint(feature_vector=[1, 2, 3])),
                FindNeighborsRequest.Query(neighbor_count=3, datapoint=IndexDatapoint(feature_vector=[4, 5, 6]))
            ],
            Parse(test_utils.read_resource("tests/unit/test_query_responses/rest_response_1.json"), FindNeighborsResponse()._pb),
            MatchResult(matches=[
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2]),
                    matching_client.MatchResult.Neighbor(id="1", distance=10.0, feature_vector=[3, 4, 5]),
                    matching_client.MatchResult.Neighbor(id="2", distance=20.0, feature_vector=[6, 7, 8]),
                ]),
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="10", distance=0.0, feature_vector=[10, 11, 12]),
                    matching_client.MatchResult.Neighbor(id="11", distance=0.1, feature_vector=[13, 14, 15]),
                    matching_client.MatchResult.Neighbor(id="12", distance=0.2, feature_vector=[16, 17, 18]),
                ])
            ])
        ),
    ])
    def test_rest_client(
            self,
            queries: List[List[float]],
            client_queries: List[FindNeighborsRequest.Query],
            client_response: FindNeighborsResponse,
            expected_response: MatchResult,
    ) -> None:
        """
        Test the REST matching client.

        :param queries: The sent into the match method.
        :param client_queries: The queries sent to the REST client (e.g. after translation)
        :param client_response: The response from the REST client (e.g. pre-parsing the response)
        :param expected_response: The expected response from the match method.

        """
        rest_client = mock(MatchServiceClient)
        when(rest_client).find_neighbors(FindNeighborsRequest(
            index_endpoint=REST_INDEX.endpoint_id,
            deployed_index_id=REST_INDEX.deployed_index_id,
            queries=client_queries,
            return_full_datapoint=True,
        )).thenReturn(client_response)

        match_client = MatchingClient.from_index(REST_INDEX)
        when(match_client)._MatchingClientREST__get_match_index_endpoint_client().thenReturn(rest_client)
        assert match_client.match(queries) == expected_response

    @parameterized.expand([
        ([], [], MatchResult()),
        (
            [[1, 2, 3]],
            [
                [
                    MatchNeighbor(id="0", distance=0.0, feature_vector=[
                                  0, 1, 2], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="1", distance=10.0, feature_vector=[
                                  3, 4, 5], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="2", distance=20.0, feature_vector=[
                                  6, 7, 8], crowding_tag="0", numeric_restricts=[]),
                ]
            ],
            MatchResult(matches=[
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2]),
                    matching_client.MatchResult.Neighbor(id="1", distance=10.0, feature_vector=[3, 4, 5]),
                    matching_client.MatchResult.Neighbor(id="2", distance=20.0, feature_vector=[6, 7, 8]),
                ])
            ])
        ),
        (
            [[1, 2, 3], [4, 5, 6]],
            [
                [
                    MatchNeighbor(id="0", distance=0.0, feature_vector=[
                                  0, 1, 2], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="1", distance=10.0, feature_vector=[
                                  3, 4, 5], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="2", distance=20.0, feature_vector=[
                                  6, 7, 8], crowding_tag="0", numeric_restricts=[]),
                ],
                [
                    MatchNeighbor(id="10", distance=0.0, feature_vector=[
                                  10, 11, 12], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="11", distance=0.1, feature_vector=[
                                  13, 14, 15], crowding_tag="0", numeric_restricts=[]),
                    MatchNeighbor(id="12", distance=0.2, feature_vector=[
                                  16, 17, 18], crowding_tag="0", numeric_restricts=[]),
                ]
            ],
            MatchResult(matches=[
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="0", distance=0.0, feature_vector=[0, 1, 2]),
                    matching_client.MatchResult.Neighbor(id="1", distance=10.0, feature_vector=[3, 4, 5]),
                    matching_client.MatchResult.Neighbor(id="2", distance=20.0, feature_vector=[6, 7, 8]),
                ]),
                matching_client.MatchResult.NearestNeighbors(neighbors=[
                    matching_client.MatchResult.Neighbor(id="10", distance=0.0, feature_vector=[10, 11, 12]),
                    matching_client.MatchResult.Neighbor(id="11", distance=0.1, feature_vector=[13, 14, 15]),
                    matching_client.MatchResult.Neighbor(id="12", distance=0.2, feature_vector=[16, 17, 18]),
                ])
            ])
        ),
    ])
    def test_grpc_client(
            self,
            queries: List[List[float]],
            client_response: List[List[MatchNeighbor]],
            expected_response: MatchResult,
    ) -> None:
        """
        Test the gRPC matching client.

        :param queries: The sent into the match method.
        :param client_response: The response from the gRPC client (e.g. pre-parsing the response)
        :param expected_response: The expected response from the match method.

        """
        grpc_client = mock(CustomMatchingEngineIndexEndpointClient)
        when(grpc_client).match(
            deployed_index_id=GRPC_INDEX.deployed_index_id,
            queries=queries,
            num_neighbors=GRPC_INDEX.num_neighbors,
        ).thenReturn(client_response)

        match_client = MatchingClient.from_index(GRPC_INDEX)
        when(match_client)._MatchingClientGRPC__get_match_index_endpoint_client().thenReturn(grpc_client)
        assert match_client.match(queries) == expected_response
