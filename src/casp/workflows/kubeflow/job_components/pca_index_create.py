from kfp import dsl

from casp.workflows.kubeflow import machine_specs


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def create_deploy_register_index(gcs_config_path: str) -> None:
    """
    Create, deploy, and register Space Vector Search index.

    :param gcs_config_path: GCS path to the index config file.
    """
    import yaml
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

            gapic_index = gca_matching_engine_index.Index(
                display_name=display_name,
                description=description,
                metadata={
                    "config": config.as_dict(),
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
        sync=False,
    )

    index.wait()
    # Deploying Index
    index_id = index.resource_name
    # Deployed Index ID cannot have dashes
    deployed_index_id = f"deployed_{config_data['display_name']}".replace("-", "_")
    endpoint = "{}-aiplatform.googleapis.com".format(config_data["location"])

    index_endpoint_client = aiplatform_v1beta1.IndexEndpointServiceClient(
        client_options={"api_endpoint": endpoint}, credentials=credentials
    )

    deploy_matching_engine_index = {
        "id": deployed_index_id,
        "display_name": deployed_index_id,
        "index": index_id,
    }

    index_endpoint_client.deploy_index(
        index_endpoint=config_data["index_endpoint_id"], deployed_index=deploy_matching_engine_index
    )

    registry_client = clients.RegistryClient()
    registry_client.register_index(
        model_name=config_data["model_name"],
        index_name=config_data["display_name"],
        deployed_index_id=deployed_index_id,
        endpoint_id=config_data["index_endpoint_id"],
        embedding_dimension=config_data["embedding_dimension"],
        num_neighbors=config_data["approximate_neighbors_count"],
    )
