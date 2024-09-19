from kfp import dsl


# @dsl.component()
def create_register_index(gcs_config_path: str) -> None:
    """
    Create, deploy, and register Space Vector Search index.

    :param gcs_config_path: GCS path to the index config file.
    """
    import typing as t

    import yaml
    from google.oauth2.service_account import Credentials
    from smart_open import open

    from casp.services import settings, utils
    from casp.workflows.kubeflow.job_components_library import clients
    from casp.workflows.kubeflow.job_components_library.matching_engine_client import CustomMatchingEngineIndex

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    credentials = Credentials.from_service_account_info(
        info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    project_id = settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")

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

    # clients.CellariumCloudInternalServiceClient.register_index(
    #     model_name=config_data["model_name"],
    #     index_name=config_data["display_name"],
    #     deployed_index_id=None,
    #     endpoint_id=config_data["index_endpoint_id"],
    #     embedding_dimension=config_data["embedding_dimension"],
    #     num_neighbors=config_data["approximate_neighbors_count"],
    # )
