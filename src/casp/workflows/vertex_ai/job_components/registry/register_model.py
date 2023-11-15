import typing as t

from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component
from kfp import dsl

from casp.workflows.vertex_ai.job_components import constants


@dsl.component(base_image=constants.DOCKER_IMAGE_NAME_CPU)
def register_embedding_model(gcs_config_path: str) -> None:
    import yaml
    from smart_open import open

    from casp.workflows.vertex_ai.job_components.services import EmbeddingModelRegistryService

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    service = EmbeddingModelRegistryService()
    service.register_embedding_model(
        model_name=config_data["model_name"],
        model_file_path=config_data["model_file_path"],
        embedding_dimension=config_data["embedding_dimension"],
    )


def create_register_embedding_model_job(gcs_config_path: str) -> t.Callable[[], t.Any]:
    register_embedding_model_job = create_custom_training_job_from_component(
        register_embedding_model,
        display_name=constants.REGISTRY_DISPLAY_NAME,
        replica_count=constants.REGISTRY_REPLICA_COUNT,
        machine_type=constants.REGISTRY_MACHINE_TYPE,
    )
    return lambda: register_embedding_model_job(gcs_config_path=gcs_config_path)
