import typing as t

from casp.clients import HTTPClient
from casp.services import settings


class RegistryClient(HTTPClient):
    BACKEND_URL: str = f"{settings.API_INTERNAL_SERVER_URL}/api"

    def register_model(
        self, model_name: str, model_file_path: str, embedding_dimension: int, bq_dataset_name: str, schema_name: str
    ) -> t.Dict[str, t.Any]:
        """
        Register a model in the Cellarium Cloud internal database using the API Internal server.

        :param model_name: Model name
        :param model_file_path: Model file path in GCS
        :param embedding_dimension: Model embedding dimension
        :param bq_dataset_name: BigQuery dataset name which was used to train the model
        :param schema_name: Schema name which was used to train the model

        :return: Dictionary with the response from the server
        """
        data = {
            "model_name": model_name,
            "model_file_path": model_file_path,
            "embedding_dimension": embedding_dimension,
            "bq_dataset_name": bq_dataset_name,
            "schema_name": schema_name,
        }
        return self.post_json(endpoint="cas-model", data=data)

    def register_index(
        self,
        model_name: str,
        index_name: str,
        num_neighbors: int,
        deployed_index_id: str,
        endpoint_id: str,
        embedding_dimension: int,
    ) -> t.Dict[str, t.Any]:
        """
        Register an index in the Cellarium Cloud internal database using the API Internal server.

        :param model_name: Model name to which the index belongs
        :param index_name: Index name
        :param num_neighbors: Number of neighbors that will be returned by the index
        :param deployed_index_id: Deployed index ID
        :param endpoint_id: Endpoint ID where the index is deployed
        :param embedding_dimension: Index embedding dimension

        :return: Dictionary with the response from the server
        """
        data = {
            "model_name": model_name,
            "index_name": index_name,
            "deployed_index_id": deployed_index_id,
            "endpoint_id": endpoint_id,
            "embedding_dimension": embedding_dimension,
            "num_neighbors": num_neighbors,
        }
        return self.post_json(endpoint="cas-index", data=data)
