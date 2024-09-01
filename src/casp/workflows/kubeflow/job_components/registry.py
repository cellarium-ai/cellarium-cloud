from kfp import dsl


# @dsl.component()
def register_embedding_model(gcs_config_path: str) -> None:
    """
    Component for registering embedding model in the registry by creating an instance in the database with the model
    information

    :param gcs_config_path: GCS path to the config file with the model information
    """
    import json
    import typing as t

    import requests
    import yaml
    from smart_open import open

    # from casp.workflows.kubeflow.job_components_library.clients import RegistryClient

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    class RegistryClient:
        def register_model(
            self,
            model_name: str,
            model_file_path: str,
            embedding_dimension: int,
            bq_dataset_name: str,
            schema_name: str,
        ):
            """
            Register a model in the Cellarium Cloud internal database using the API Internal server.

            :param model_name: Model name
            :param model_file_path: Model file path in GCS
            :param embedding_dimension: Model embedding dimension
            :param bq_dataset_name: BigQuery dataset name which was used to train the model
            :param schema_name: Schema name which was used to train the model
            """
            url = "https://cas-api-internal-1-4-2-dev-vi7nxpvk7a-uc.a.run.app/api/cas-model"

            data = {
                "model_name": model_name,
                "model_file_path": model_file_path,
                "embedding_dimension": embedding_dimension,
                "admin_use_only": True,
                "bq_dataset_name": bq_dataset_name,
                "schema_name": schema_name,
            }
            headers = {"Content-Type": "application/json"}

            response = requests.post(url=url, data=json.dumps(data), headers=headers)
            # Print the response
            print(f"Status Code: {response.status_code}")
            print(f"Response Text: {response.text}")

        def register_index(
            self,
            model_name: str,
            index_name: str,
            num_neighbors: int,
            endpoint_id: str,
            embedding_dimension: int,
            deployed_index_id: t.Optional[str] = None,
        ):
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
            url = "https://cas-api-internal-1-4-2-dev-vi7nxpvk7a-uc.a.run.app/api/cas-index"

            data = {
                "model_name": model_name,
                "index_name": index_name,
                "embedding_dimension": embedding_dimension,
                "endpoint_id": endpoint_id,
                "deployed_index_id": deployed_index_id,
                "num_neighbors": num_neighbors,
                "is_grpc": True,
            }
            headers = {"Content-Type": "application/json"}
            # Make the POST request with Basic Auth
            response = requests.post(url, data=json.dumps(data), headers=headers)

            # Print the response
            print(response.status_code)
            print(response.text)

    client = RegistryClient()
    client.register_model(
        model_name=config_data["model_name"],
        model_file_path=config_data["model_file_path"],
        embedding_dimension=config_data["embedding_dimension"],
        bq_dataset_name=config_data["bq_dataset_name"],
        schema_name=config_data["schema_name"],
    )
