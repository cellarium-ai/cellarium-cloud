import typing as t

from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component
from kfp import dsl

from casp.workflows.vertex_ai.job_components import constants


@dsl.component(base_image=constants.DOCKER_IMAGE_NAME_CUDA)
def cellarium_pca_train_component(gcs_config_path: str):
    """
    Cellarium ML CLI component for running incremental PCA training.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK") is not None:
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["incremental_pca", "fit", "--config", gcs_config_path])


def create_pca_train_job(gcs_config_path: str) -> t.Callable[[], t.Any]:
    """
    Create a custom training Google Vertex AI job for running incremental PCA training.

    :param gcs_config_path: Config file path on GCS.

    :return: Callable custom training job.
    """
    pca_train_job = create_custom_training_job_from_component(
        cellarium_pca_train_component,
        display_name=constants.TRAIN_DISPLAY_NAME,
        replica_count=constants.TRAIN_REPLICA_COUNT,
        machine_type=constants.TRAIN_MACHINE_TYPE,
        accelerator_type=constants.TRAIN_ACCELERATOR_TYPE,
        accelerator_count=constants.TRAIN_ACCELERATOR_COUNT,
    )
    return lambda: pca_train_job(gcs_config_path=gcs_config_path)
