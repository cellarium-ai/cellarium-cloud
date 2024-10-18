import typer

from casp.workflows.kubeflow.pipelines import pipelines
from casp.workflows.kubeflow import kubeflow_command_registry

typer_app = typer.Typer()

# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="calculate-metrics",
#     pipeline_func=pipelines.calculate_metrics_pipeline,
#     display_name="Calculate Metrics",
#     help_text="Run a pipeline to calculate metrics",
# )

# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="generate-cas-outputs",
#     pipeline_func=pipelines.generate_cas_outputs_pipeline,
#     display_name="Generate CAS Outputs",
#     help_text="Run a pipeline to generate CAS outputs",
# )


kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="bq-ops-prepare-and-extract",
    pipeline_func=pipelines.bq_ops_prepare_and_extract,
    display_name="BQ Ops prepare and extract",
    help_text="BQ Ops prepare extract tables and process extract",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="bq-ops-prepare-extract",
    pipeline_func=pipelines.bq_ops_prepare_extract,
    display_name="BQ Ops prepare extract",
    help_text="BQ Ops prepare extract tables",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="bq-ops-extract",
    pipeline_func=pipelines.bq_ops_extract,
    display_name="BQ Ops extract",
    help_text="BQ Ops extract data",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="summary-stats-train",
    pipeline_func=pipelines.summary_stats_train_pipeline,
    display_name="Train Summary Stats",
    help_text="Train Summary Statis",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="pca-train",
    pipeline_func=pipelines.pca_train_pipeline,
    display_name="Train PCA Model",
    help_text="Train PCA Model",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="pca-resize-embed",
    pipeline_func=pipelines.pca_resize_and_embed_pipeline,
    display_name="Resize PCA model and embed data",
    help_text="Resize PCA model and embed data",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="pca-resize-full",
    pipeline_func=pipelines.pca_resize_full_cycle_pipeline,
    display_name="Resize PCA model and create all resources",
    help_text="Resize PCA model and create all resources",
)


kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="benchmarking-generate-outputs",
    pipeline_func=pipelines.generate_cas_outputs_pipeline,
    display_name="Generate CAS Outputs for Benchmarking",
    help_text="Generate CAS Outputs for Benchmarking",
)

kubeflow_command_registry.register_pipeline_as_command(
    typer_app=typer_app,
    name="benchmarking-calculate-metrics",
    pipeline_func=pipelines.calculate_metrics_pipeline,
    display_name="benchmarking_calculate_metrics",
    help_text="Benchmark CAS results",
)

# @typer_app.command()
# def summary_stats_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run summary stats train pipeline in parallel.
#     This pipeline will train onepass_mean_var_std and tdigest models and save them in a GCS bucket.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.summary_stats_train_pipeline,
#         pipeline_display_name="summary_stats_train_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )


# @typer_app.command()
# def bq_ops_prepare_extract(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_prepare_extract,
#         pipeline_display_name="Prepare extract table in BigQuery",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def bq_ops_extract(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_extract,
#         pipeline_display_name="Extract files from BigQuery",
#         pipeline_kwargs=config_paths,
#     )


# @typer_app.command()
# def pca_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA train pipeline in parallel.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipeline_func=pipelines.train_pipeline,
#         pipeline_display_name="ipca_train",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )


#
#
# @typer_app.command()
# def pca_embed(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA embed pipeline in parallel.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipeline_func=pipelines.embed_pipeline,
#         pipeline_display_name="pca_embed",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )


#
# @typer_app.command()
# def pca_train_embed(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA train and embed pipelines in sequence for multiple models simultaneously.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.train_embed_pipeline,
#         pipeline_display_name="ipca_train_embed_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def pca_deploy_vector_search_index(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA deploy vector search index pipeline in parallel.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.pca_deploy_index_pipeline,
#         pipeline_display_name="pca_deploy_vector_search_index_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def pca_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA full cycle pipeline in parallel.
#
#     Full cycle pipeline consists of the following steps:
#     - Train PCA model
#     - Embed data using the trained PCA model
#     - Register the trained PCA model in Cellarium Cloud admin
#     - Create Vector search index using the embedded data
#     - Deploy the vector search index to the endpoint
#     - Register the vector search index in Cellarium Cloud admin
#
#     .. note::
#         Trained models will be registered with `admin_use_only=True` flag.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.pca_full_cycle_pipeline,
#         pipeline_display_name="pca_full_cycle_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def pca_resize_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run PCA resize and save pipeline in parallel.
#
#     Full cycle resize pipeline consists of the following steps:
#     - Take a trained PCA model as a base, slice its V_kg and S_k to the new embedding dimension and save the new model.
#     - Embed data using a new PCA model
#     - Register the resized PCA model in Cellarium Cloud admin
#     - Create Vector search index using the embedded data
#     - Deploy the vector search index to the endpoint
#     - Register the vector search index in Cellarium Cloud admin
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.pca_resize_full_cycle_pipeline,
#         pipeline_display_name="pca_resize_full_cycle_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def summary_stats_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run summary stats train pipeline in parallel.
#     This pipeline will train onepass_mean_var_std and tdigest models and save them in a GCS bucket.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.summary_stats_train_pipeline,
#         pipeline_display_name="summary_stats_train_parallel",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def logistic_regression_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run logistic regression train pipeline in parallel.
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     pipeline_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.logistic_regression_train_pipeline,
#         pipeline_display_name="logistic_regression_train_parallel",
#         pipeline_kwargs=pipeline_kwargs,
#     )
#
#
# @typer_app.command()
# def pca_train_base_resize_full_cycle(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.pca_train_base_and_resize_full_cycle_pipeline,
#         pipeline_display_name="pca_train_base_resize_full_cycle_parallel",
#         pipeline_kwargs=pipeline_as_components_kwargs,
#     )
#
#
# @typer_app.command()
# def summary_stats_and_logistic_regression_train(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.summary_stats_lr_train_pipeline,
#         pipeline_display_name="summary_stats_and_logistic_regression_train_parallel",
#         pipeline_kwargs=pipeline_as_components_kwargs,
#     )
#
#
# @typer_app.command()
# def generate_cas_outputs(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run Generate CAS Outputs
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.generate_cas_outputs_pipeline,
#         pipeline_display_name="generate_cas_outputs",
#         pipeline_kwargs=pipeline_as_components_kwargs,
#     )
#
#
# @typer_app.command()
# def calculate_metrics(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     """
#     Run Calculate Metrics
#
#     :param config_yaml_path: Path to the local YAML config file containing a list of configs for each run.
#     """
#     pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.calculate_metrics_pipeline,
#         pipeline_display_name="calculate_metrics",
#         pipeline_kwargs=pipeline_as_components_kwargs,
#     )
#
#
# @typer_app.command()
# def bq_ops_create_avro_files(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_create_avro_files,
#         pipeline_display_name="Create avro files for BigQuery ingest",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def bq_ops_ingest_data(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_ingest_data,
#         pipeline_display_name="Ingest data to BigQuery",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def bq_ops_create_avro_files_and_ingest(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     pipeline_as_components_kwargs = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_create_avro_files_and_ingest_data,
#         pipeline_display_name="bq_ops_create_avro_files_and_ingest",
#         pipeline_kwargs=pipeline_as_components_kwargs,
#     )
#
#
# @typer_app.command()
# def bq_ops_precalculate_fields(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_precalculate_fields,
#         pipeline_display_name="Precalculate fields in BigQuery",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def bq_ops_prepare_extract(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_prepare_extract,
#         pipeline_display_name="Prepare extract table in BigQuery",
#         pipeline_kwargs={"pipeline_config_paths": config_paths},
#     )
#
#
# @typer_app.command()
# def bq_ops_extract(config_yaml_path: Annotated[str, typer.Option()]) -> None:
#     config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
#     submit_pipeline(
#         pipelines.bq_ops_extract,
#         pipeline_display_name="Extract files from BigQuery",
#         pipeline_kwargs=config_paths,
#     )


if __name__ == "__main__":
    typer_app()
