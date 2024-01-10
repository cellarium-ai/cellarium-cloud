from kfp import dsl

from casp.workflows.kubeflow import machine_specs


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CUDA)
def train(gcs_config_path: str):
    """
    Cellarium ML CLI component for running incremental PCA training.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK") is not None:
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["incremental_pca", "fit", "--config", gcs_config_path])


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CUDA)
def embed(gcs_config_path: str) -> None:
    """
    Cellarium ML CLI component for running incremental PCA Data embedding.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK"):
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["incremental_pca", "predict", "--config", gcs_config_path])
