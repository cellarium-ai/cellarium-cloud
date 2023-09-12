import typing as t

import grpc
from google.cloud import aiplatform
from google.cloud.aiplatform.matching_engine._protos import match_service_pb2, match_service_pb2_grpc
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor

if t.TYPE_CHECKING:
    from google.auth import credentials as auth_credentials
    from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import Namespace

GRPC_MESSAGE_MAX_SIZE_DEFAULT = 50000000  # 5mb


class MatchingEngineIndexEndpointGRPCOptionsExposed(aiplatform.MatchingEngineIndexEndpoint):
    """
    Matching Engine index endpoint resource for Vertex AI. Inherited from `aiplatform.MatchingEngineIndexEndpoint`.
    This class provides additional flexibility to increase gRPC communication protocol message maximum size.
    """

    def __init__(
        self,
        index_endpoint_name: str,
        project: t.Optional[str] = None,
        location: t.Optional[str] = None,
        credentials: t.Optional["auth_credentials.Credentials"] = None,
        grpc_message_max_size: int = GRPC_MESSAGE_MAX_SIZE_DEFAULT,
    ):
        """
        Retrieves an existing index endpoint given a name or ID.

        Example Usage:

           my_index_endpoint = aiplatform.MatchingEngineIndexEndpoint(
               index_endpoint_name='projects/123/locations/us-central1/index_endpoint/my_index_id',
               grpc_message_max_size=10000
           )
           or
           my_index_endpoint = aiplatform.MatchingEngineIndexEndpoint(
               index_endpoint_name='my_index_endpoint_id'
           )

        :param index_endpoint_name: A fully-qualified index endpoint resource name or a index ID.
           Example: "projects/123/locations/us-central1/index_endpoints/my_index_id"
           or "my_index_id" when project and location are initialized or passed.
        :param project:
           Project to retrieve index endpoint from. If not set, project
           set in `aiplatform.init` will be used.
        :param location:
           Location to retrieve index endpoint from. If not set, location
           set in `aiplatform.init` will be used.
        :param credentials:
           Custom credentials to use to retrieve this IndexEndpoint. Overrides
           credentials set in aiplatform.init.
        :param grpc_message_max_size: gRPC message maximum size
        """
        self.grpc_message_max_size = grpc_message_max_size
        super().__init__(
            index_endpoint_name=index_endpoint_name, project=project, location=location, credentials=credentials
        )

    def match(
        self,
        deployed_index_id: str,
        queries: t.List[t.List[float]],
        num_neighbors: int = 1,
        filter: t.Optional[t.List["Namespace"]] = None,
    ) -> t.List[t.List["MatchNeighbor"]]:
        """
        Method is overriden to be capable of larger request sizes. For more info refer to:
        `aiplatform.MatchingEngineIndexEndpoint` class.

        :param deployed_index_id: The ID of the DeployedIndex to match the queries against.
        :param queries: A list of queries. Each query is a list of floats, representing a single embedding.
        :param num_neighbors: The number of nearest neighbors to be retrieved from database for each query.
        :param filter: A list of Namespaces for filtering the matching results.
                For example, [Namespace("color", ["red"], []), Namespace("shape", [], ["squared"])] will match
                datapoints that satisfy "red color" but not include datapoints with "squared shape".
                Please refer to https://cloud.google.com/vertex-ai/docs/matching-engine/filtering#json for more detail.
        :return:  A list of nearest neighbors for each query.
        """
        # ==== Overriden fix: make an input argument immutable ====
        if filter is None:
            filter = []
        # ==== overriden fix ends ====
        # Find the deployed index by id
        deployed_indexes = [
            deployed_index for deployed_index in self.deployed_indexes if deployed_index.id == deployed_index_id
        ]

        if not deployed_indexes:
            raise RuntimeError(f"No deployed index with id '{deployed_index_id}' found")

        # Retrieve server ip from deployed index
        server_ip = deployed_indexes[0].private_endpoints.match_grpc_address

        # Set up channel and stub
        # ==== Here is overridden part ====
        channel = grpc.insecure_channel(
            "{}:10000".format(server_ip),
            options=[
                ("grpc.max_send_message_length", self.grpc_message_max_size),
                ("grpc.max_receive_message_length", self.grpc_message_max_size),
            ],
        )
        # ==== overridden part ends ====
        stub = match_service_pb2_grpc.MatchServiceStub(channel)

        # Create the batch match request
        batch_request = match_service_pb2.BatchMatchRequest()
        batch_request_for_index = match_service_pb2.BatchMatchRequest.BatchMatchRequestPerIndex()
        batch_request_for_index.deployed_index_id = deployed_index_id
        requests = []
        for query in queries:
            request = match_service_pb2.MatchRequest(
                num_neighbors=num_neighbors,
                deployed_index_id=deployed_index_id,
                float_val=query,
            )
            for namespace in filter:
                restrict = match_service_pb2.Namespace()
                restrict.name = namespace.name
                restrict.allow_tokens.extend(namespace.allow_tokens)
                restrict.deny_tokens.extend(namespace.deny_tokens)
                request.restricts.append(restrict)
            requests.append(request)

        batch_request_for_index.requests.extend(requests)
        batch_request.requests.append(batch_request_for_index)

        # Perform the request
        response = stub.BatchMatch(batch_request)

        # Wrap the results in MatchNeighbor objects and return
        return [
            [MatchNeighbor(id=neighbor.id, distance=neighbor.distance) for neighbor in embedding_neighbors.neighbor]
            for embedding_neighbors in response.responses[0].responses
        ]
