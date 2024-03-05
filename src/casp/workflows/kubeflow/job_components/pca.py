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


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def resize_and_save(gcs_config_path: str) -> None:
    import torch
    import yaml
    from smart_open import open

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    with open(config_data["base_model_path"], "rb") as f:
        model_ckpt = torch.load(f, map_location="cpu")

    new_embedding_dimension = config_data["new_embedding_dimension"]

    state_dict_keys = model_ckpt["state_dict"].keys()
    model_pipeline_number = 0
    if "pipeline.1.V_kg" in state_dict_keys:
        model_pipeline_number = 1
    elif "pipeline.2.V_kg" in state_dict_keys:
        model_pipeline_number = 2
    elif "pipeline.3.V_kg" in state_dict_keys:
        model_pipeline_number = 3

    model_ckpt["state_dict"][f"pipeline.{model_pipeline_number}.V_kg"] = model_ckpt["state_dict"][
        f"pipeline.{model_pipeline_number}.V_kg"
    ][:new_embedding_dimension, :].clone()
    model_ckpt["state_dict"][f"pipeline.{model_pipeline_number}.S_k"] = model_ckpt["state_dict"][
        f"pipeline.{model_pipeline_number}.S_k"
    ][:new_embedding_dimension].clone()

    model_ckpt["hyper_parameters"]["model"]["init_args"]["n_components"] = new_embedding_dimension

    with open(config_data["new_model_path"], "wb") as file:
        torch.save(model_ckpt, file)
