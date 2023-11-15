import json

import typer
from typing_extensions import Annotated

from casp.workflows.vertex_ai.pipelines import pca_pipelines
from casp.workflows.vertex_ai.pipelines.config_management import create_pca_configs_from_yaml
from casp.workflows.vertex_ai.pipelines.utils import submit_cellarium_ml_pipeline

typer_app = typer.Typer()


@typer_app.command()
def train(gcs_config_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA train pipeline.

    :param gcs_config_path: GCS path to the PCA train config file.
    """
    submit_cellarium_ml_pipeline(
        pipeline_func=pca_pipelines.train_pipeline,
        pipeline_display_name="ipca_train",
        pipeline_kwargs={"gcs_config_path": gcs_config_path},
    )


@typer_app.command()
def embed(gcs_config_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA embed pipeline.

    :param gcs_config_path: GCS path to the PCA embed config file.
    """
    submit_cellarium_ml_pipeline(
        pipeline_func=pca_pipelines.embed_pipeline,
        pipeline_display_name="ipca_embed",
        pipeline_kwargs={"gcs_config_path": gcs_config_path},
    )


@typer_app.command()
def train_embed(
    train_gcs_config_path: Annotated[str, typer.Option()], embed_gcs_config_path: Annotated[str, typer.Option()]
) -> None:
    """
    Run PCA train and embed pipelines in sequence.

    :param train_gcs_config_path: GCS path to the PCA train config file.
    :param embed_gcs_config_path: GCS path to the PCA embed config file.
    """
    pipeline_kwargs = {"train_gcs_config_path": train_gcs_config_path, "embed_gcs_config_path": embed_gcs_config_path}
    submit_cellarium_ml_pipeline(
        pipeline_func=pca_pipelines.train_embed_pipeline,
        pipeline_display_name="ipca_train_embed",
        pipeline_kwargs=pipeline_kwargs,
    )


@typer_app.command()
def train_embed_bulk(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    config_paths = create_pca_configs_from_yaml(config_yaml_path)

    # DSL Pipelines require JSON strings as input
    config_paths = [json.dumps(config) for config in config_paths]

    submit_cellarium_ml_pipeline(
        pca_pipelines.train_embed_bulk_pipeline,
        pipeline_display_name="ipca_train_embed_bulk",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


if __name__ == "__main__":
    typer_app()
