from typing import Dict, Optional, Sequence, Tuple

from google.auth import credentials as auth_credentials
from google.cloud import aiplatform
from google.cloud.aiplatform import base, initializer, utils
from google.cloud.aiplatform.compat.types import encryption_spec as gca_encryption_spec
from google.cloud.aiplatform.compat.types import matching_engine_index as gca_matching_engine_index
from google.cloud.aiplatform.matching_engine import matching_engine_index_config
from google.cloud.aiplatform_v1 import IndexEndpointServiceClient
from google.cloud.aiplatform_v1.types import DeployedIndex, DeployIndexRequest, UndeployIndexRequest

from casp.services import utils

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
    ) -> "aiplatform.MatchingEngineIndex":
        """
        Override creation from :class:`aiplatform.MatchingEngineIndex`. Allow using featureNormType in
        `request_metadata` as parent class doesn't handel this parameter.
        """
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

    @staticmethod
    def deploy_index(region, project_id, index_id, index_endpoint_id, display_name, deployed_index_id):
        # Initialize Vertex AI SDK
        aiplatform.init(
            project=project_id,
            location=region,
        )

        # Define variables
        index_endpoint_name = f"projects/{project_id}/locations/{region}/indexEndpoints/{index_endpoint_id}"
        index_name = f"projects/{project_id}/locations/{region}/indexes/{index_id}"

        # Create client
        client = IndexEndpointServiceClient(client_options={"api_endpoint": f"{region}-aiplatform.googleapis.com"})

        # Create deployed index
        deployed_index = DeployedIndex(
            id=deployed_index_id,
            index=index_name,
            display_name=display_name,
            automatic_resources={"min_replica_count": 1, "max_replica_count": 4},
        )

        # Prepare deployment request
        deploy_index_request = DeployIndexRequest(
            index_endpoint=index_endpoint_name,
            deployed_index=deployed_index,
        )

        print("Starting index deployment...")
        operation = client.deploy_index(request=deploy_index_request)

        # Wait for the operation to complete
        print("Waiting for the deployment to complete...")
        operation.result(timeout=7200)  # Timeout after 2 hours

        # Check if the operation was successful
        if operation.done():
            print("Index deployment completed successfully.")
        else:
            print("Index deployment failed. Please check the Google Cloud Console for details.")

    @staticmethod
    def undeploy_index(region, project_id, index_endpoint_id, deployed_index_id):
        # Create client
        client = IndexEndpointServiceClient(client_options={"api_endpoint": f"{region}-aiplatform.googleapis.com"})
        index_endpoint_name = f"projects/{project_id}/locations/{region}/indexEndpoints/{index_endpoint_id}"
        # Prepare undeployment request
        undeploy_index_request = UndeployIndexRequest(
            index_endpoint=index_endpoint_name,
            deployed_index_id=deployed_index_id,
        )

        # Undeploy the index for the operation to complete
        print("Starting index undeployment...")
        _ = client.undeploy_index(request=undeploy_index_request)
