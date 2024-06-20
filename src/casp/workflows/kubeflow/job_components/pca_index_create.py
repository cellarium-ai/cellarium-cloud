from kfp import dsl


@dsl.component()
def create_register_index(gcs_config_path: str) -> None:
    """
    Create, deploy, and register Space Vector Search index.

    :param gcs_config_path: GCS path to the index config file.
    """
    import yaml
    import json
    import typing as t
    import requests
    from google.cloud import aiplatform, aiplatform_v1beta1
    from smart_open import open

    from casp.services import utils
    from casp.workflows.kubeflow.job_components_library import clients

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    credentials, project_id = utils.get_google_service_credentials()
    from typing import Dict, List, Optional, Sequence, Tuple

    from google.auth import credentials as auth_credentials
    from google.cloud.aiplatform import base, initializer, utils
    from google.cloud.aiplatform.compat.types import encryption_spec as gca_encryption_spec
    from google.cloud.aiplatform.compat.types import index_service as gca_index_service
    from google.cloud.aiplatform.compat.types import (
        matching_engine_deployed_index_ref as gca_matching_engine_deployed_index_ref,
    )
    from google.cloud.aiplatform.compat.types import matching_engine_index as gca_matching_engine_index
    from google.cloud.aiplatform.matching_engine import matching_engine_index_config

    _LOGGER = base.Logger(__name__)
    _INDEX_UPDATE_METHOD_TO_ENUM_VALUE = {
        "STREAM_UPDATE": gca_matching_engine_index.Index.IndexUpdateMethod.STREAM_UPDATE,
        "BATCH_UPDATE": gca_matching_engine_index.Index.IndexUpdateMethod.BATCH_UPDATE,
    }

    class CustomMatchingEngineIndex(aiplatform.MatchingEngineIndex):
        @classmethod
        @base.optional_sync()
        def _create(
            cls,
            display_name: str,
            contents_delta_uri: str,
            config: matching_engine_index_config.MatchingEngineIndexConfig,
            description: Optional[str] = None,
            labels: Optional[Dict[str, str]] = None,
            project: Optional[str] = None,
            location: Optional[str] = None,
            credentials: Optional[auth_credentials.Credentials] = None,
            request_metadata: Optional[Sequence[Tuple[str, str]]] = (),
            sync: bool = True,
            index_update_method: Optional[str] = None,
            encryption_spec_key_name: Optional[str] = None,
        ) -> "MatchingEngineIndex":

            index_update_method_enum = None
            if index_update_method in _INDEX_UPDATE_METHOD_TO_ENUM_VALUE:
                index_update_method_enum = _INDEX_UPDATE_METHOD_TO_ENUM_VALUE[index_update_method]

            config = config.as_dict()
            request_metadata_keys = [k for k, _ in request_metadata]

            if "featureNormType" in request_metadata_keys:
                feature_norm_type = [x[-1] for x in request_metadata if x[0] == "featureNormType"][0]
                # pop the featureNormType from request_metadata
                request_metadata = [x for x in request_metadata if x[0] != "featureNormType"]
                config["featureNormType"] = feature_norm_type

            gapic_index = gca_matching_engine_index.Index(
                display_name=display_name,
                description=description,
                metadata={
                    "config": config,
                    "contentsDeltaUri": contents_delta_uri,
                },
                index_update_method=index_update_method_enum,
            )

            if encryption_spec_key_name:
                encryption_spec = gca_encryption_spec.EncryptionSpec(kms_key_name=encryption_spec_key_name)
                gapic_index.encryption_spec = encryption_spec

            if labels:
                utils.validate_labels(labels)
                gapic_index.labels = labels

            api_client = cls._instantiate_client(location=location, credentials=credentials)

            create_lro = api_client.create_index(
                parent=initializer.global_config.common_location_path(project=project, location=location),
                index=gapic_index,
                metadata=request_metadata,
            )

            _LOGGER.log_create_with_lro(cls, create_lro)

            created_index = create_lro.result(timeout=10000)

            _LOGGER.log_create_complete(cls, created_index, "index")

            index_obj = cls(
                index_name=created_index.name,
                project=project,
                location=location,
                credentials=credentials,
            )

            return index_obj

    request_metadata = {"featureNormType": config_data["feature_norm_type"]}
    request_metadata = [(k, v) for k, v in request_metadata.items()]
    # Creating Index
    index = CustomMatchingEngineIndex.create_tree_ah_index(
        display_name=config_data["display_name"],
        contents_delta_uri=config_data["contents_delta_uri"],
        dimensions=config_data["embedding_dimension"],
        approximate_neighbors_count=config_data["approximate_neighbors_count"],
        project=project_id,
        location=config_data["location"],
        credentials=credentials,
        distance_measure_type=config_data["distance_measure_type"],
        leaf_node_embedding_count=config_data["leaf_node_embedding_count"],
        request_metadata=request_metadata,
        sync=False,
    )

    index.wait()

    #### Deploying Index ###
    # index_id = index.resource_name
    # # Deployed Index ID cannot have dashes
    # deployed_index_id = f"deployed_{config_data['display_name']}".replace("-", "_")
    # endpoint = "{}-aiplatform.googleapis.com".format(config_data["location"])
    #
    # index_endpoint_client = aiplatform_v1beta1.IndexEndpointServiceClient(
    #     client_options={"api_endpoint": endpoint}, credentials=credentials
    # )
    #
    # deploy_matching_engine_index = {
    #     "id": deployed_index_id,
    #     "display_name": deployed_index_id,
    #     "index": index_id,
    # }
    #
    # index_endpoint_client.deploy_index(
    #     index_endpoint=config_data["index_endpoint_id"], deployed_index=deploy_matching_engine_index
    # )
    ########################

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

    deployed_index_id = None
    registry_client = RegistryClient()
    registry_client.register_index(
        model_name=config_data["model_name"],
        index_name=config_data["display_name"],
        deployed_index_id=deployed_index_id,
        endpoint_id=config_data["index_endpoint_id"],
        embedding_dimension=config_data["embedding_dimension"],
        num_neighbors=config_data["approximate_neighbors_count"],
    )
