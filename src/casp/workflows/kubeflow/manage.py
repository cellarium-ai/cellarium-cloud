import typer
from typing_extensions import Annotated

from casp.workflows.kubeflow import config_management
from casp.workflows.kubeflow.pipelines import pipelines
from casp.workflows.kubeflow.pipelines.utils import submit_pipeline

typer_app = typer.Typer()


@typer_app.command()
def pca_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA train pipeline in parallel.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipeline_func=pipelines.train_pipeline,
        pipeline_display_name="ipca_train",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_embed(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA embed pipeline in parallel.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipeline_func=pipelines.embed_pipeline,
        pipeline_display_name="pca_embed",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_train_embed(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA train and embed pipelines in sequence for multiple models simultaneously.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.train_embed_pipeline,
        pipeline_display_name="ipca_train_embed_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_deploy_vector_search_index(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA deploy vector search index pipeline in parallel.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.pca_deploy_index_pipeline,
        pipeline_display_name="pca_deploy_vector_search_index_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA full cycle pipeline in parallel.

    Full cycle pipeline consists of the following steps:
    - Train PCA model
    - Embed data using the trained PCA model
    - Register the trained PCA model in Cellarium Cloud admin
    - Create Vector search index using the embedded data
    - Deploy the vector search index to the endpoint
    - Register the vector search index in Cellarium Cloud admin

    .. note::
        Trained models will be registered with `admin_use_only=True` flag.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.pca_full_cycle_pipeline,
        pipeline_display_name="pca_full_cycle_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_resize_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run PCA resize and save pipeline in parallel.

    Full cycle resize pipeline consists of the following steps:
    - Take a trained PCA model as a base, slice its V_kg and S_k to the new embedding dimension and save the new model.
    - Embed data using a new PCA model
    - Register the resized PCA model in Cellarium Cloud admin
    - Create Vector search index using the embedded data
    - Deploy the vector search index to the endpoint
    - Register the vector search index in Cellarium Cloud admin

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.pca_resize_full_cycle_pipeline,
        pipeline_display_name="pca_resize_full_cycle_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def summary_stats_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run summary stats train pipeline in parallel.
    This pipeline will train onepass_mean_var_std and tdigest models and save them in a GCS bucket.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.summary_stats_train_pipeline,
        pipeline_display_name="summary_stats_train_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def logistic_regression_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run logistic regression train pipeline in parallel.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.logistic_regression_train_pipeline,
        pipeline_display_name="logistic_regression_train_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


@typer_app.command()
def pca_train_base_resize_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
    submit_pipeline(
        pipelines.pca_train_base_and_resize_full_cycle_pipeline,
        pipeline_display_name="pca_train_base_resize_full_cycle_parallel",
        pipeline_kwargs=pipeline_as_components_kwargs,
    )


@typer_app.command()
def summary_stats_and_logistic_regression_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
    submit_pipeline(
        pipelines.summary_stats_lr_train_pipeline,
        pipeline_display_name="summary_stats_and_logistic_regression_train_parallel",
        pipeline_kwargs=pipeline_as_components_kwargs,
    )


@typer_app.command()
def benchmark_cas(config_yaml_path: Annotated[str, typer.Option()]) -> None:
    """
    Run the benchmarking for CAS.

    :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
    """
    config_paths = config_management.create_configs(config_yaml_path)
    submit_pipeline(
        pipelines.benchmark_cas_pipeline,
        pipeline_display_name="benchmark_cas_parallel",
        pipeline_kwargs={"pipeline_config_paths": config_paths},
    )


if __name__ == "__main__":
    typer_app()
