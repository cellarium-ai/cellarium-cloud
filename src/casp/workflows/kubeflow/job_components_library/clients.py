import typing as t

import requests
from requests.auth import HTTPBasicAuth


class RegistryClient:

    def register_model(
        self,
        model_name: str,
        model_file_path: str,
        embedding_dimension: int,
        bq_dataset_name: str,
        schema_name: str,
        basic_auth_username: str,
        basic_auth_password: str,
    ):
        """
        Register a model in the Cellarium Cloud internal database using the API Internal server.

        :param model_name: Model name
        :param model_file_path: Model file path in GCS
        :param embedding_dimension: Model embedding dimension
        :param bq_dataset_name: BigQuery dataset name which was used to train the model
        :param schema_name: Schema name which was used to train the model
        """
        url = "https://cellarium-june-release-cas-admin-vi7nxpvk7a-uc.a.run.app/casmodel/new/?url=/casmodel/"

        form_data = {
            "model_name": model_name,
            "model_file_path": model_file_path,
            "embedding_dimension": embedding_dimension,
            "admin_use_only": "y",
            "bq_dataset_name": bq_dataset_name,
            "schema_name": schema_name,
        }

        headers = {"Content-Type": "multipart/form-data"}
        basic_auth = HTTPBasicAuth(username=basic_auth_username, password=basic_auth_password)
        response = requests.post(url=url, files=form_data, headers=headers, auth=basic_auth)

        print(response.status_code)
        print(response.text)

    def register_index(
        self,
        model_name: str,
        index_name: str,
        num_neighbors: int,
        deployed_index_id: str,
        endpoint_id: str,
        embedding_dimension: int,
        basic_auth_username: str,
        basic_auth_password: str,
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
        url = (
            "https://cellarium-june-release-cas-admin-vi7nxpvk7a-uc.a.run.app"
            "/casmatchingengineindex/new/?url=/casmatchingengineindex/"
        )

        form_data = {
            "model": "70",
            "index_name": index_name,
            "embedding_dimension": str(embedding_dimension),
            "endpoint_id": endpoint_id,
            "deployed_index_id": deployed_index_id,
            "num_neighbors": str(num_neighbors),
            "is_grpc": "y",
        }

        # Basic authentication credentials
        basic_auth = HTTPBasicAuth(username=basic_auth_username, password=basic_auth_password)
        # Make the POST request with Basic Auth
        response = requests.post(url, files=form_data, auth=basic_auth)

        # Print the response
        print(response.status_code)
        print(response.text)
