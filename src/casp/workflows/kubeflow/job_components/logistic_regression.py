from kfp import dsl

from casp.workflows.kubeflow import machine_specs


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CUDA)
def train_component(gcs_config_path: str):
    """
    Cellarium ML CLI component for running Logistic Regression model training.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK") is not None:
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["logistic_regression", "fit", "--config", gcs_config_path])
