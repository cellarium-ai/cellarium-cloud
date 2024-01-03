import typing as t

from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component
from kfp import dsl

from casp.workflows.vertex_ai.job_components import constants


@dsl.component(base_image=constants.DOCKER_IMAGE_NAME_CPU)
def create_deploy_register_index(gcs_config_path: str) -> None:
    """
    Create, deploy, and register Space Vector Search index.

    :param gcs_config_path: GCS path to the index config file.
    """
    import yaml
    from google.cloud import aiplatform, aiplatform_v1beta1
    from smart_open import open

    from casp.services import utils
    from casp.workflows.vertex_ai.job_components import clients

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    credentials, project_id = utils.get_google_service_credentials()

    # Creating Index
    index = aiplatform.MatchingEngineIndex.create_tree_ah_index(
        display_name=config_data["display_name"],
        contents_delta_uri=config_data["contents_delta_uri"],
        dimensions=config_data["embedding_dimension"],
        approximate_neighbors_count=config_data["approximate_neighbors_count"],
        project=project_id,
        location=config_data["location"],
        credentials=credentials,
        distance_measure_type=config_data["distance_measure_type"],
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


def create_index_create_job(gcs_config_path: str) -> t.Callable[[], t.Any]:
    """
    Create a custom Google Vertex AI job for creating, deploying, and registering Space Vector Search Index

    :param gcs_config_path: GCS path to the config file with the model information

    :return: Callable custom training job.
    """
    create_deploy_register_index_job = create_custom_training_job_from_component(
        create_deploy_register_index,
        display_name=constants.INDEX_CREATE_DISPLAY_NAME,
    )
    return lambda: create_deploy_register_index_job(gcs_config_path=gcs_config_path)
