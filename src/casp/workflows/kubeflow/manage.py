# import typer
import click
from casp.workflows.kubeflow.pipelines import pipelines
from casp.workflows.kubeflow import kubeflow_command_registry

# typer_app = typer.Typer()


#
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="bq-ops-prepare-and-extract",
#     pipeline_func=pipelines.bq_ops_prepare_and_extract,
#     display_name="BQ Ops prepare and extract",
#     help_text="BQ Ops prepare extract tables and process extract",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="bq-ops-prepare-extract",
#     pipeline_func=pipelines.bq_ops_prepare_extract,
#     display_name="BQ Ops prepare extract",
#     help_text="BQ Ops prepare extract tables",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="bq-ops-extract",
#     pipeline_func=pipelines.bq_ops_extract,
#     display_name="BQ Ops extract",
#     help_text="BQ Ops extract data",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="summary-stats-train",
#     pipeline_func=pipelines.summary_stats_train_pipeline,
#     display_name="Train Summary Stats",
#     help_text="Train Summary Statis",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="pca-train",
#     pipeline_func=pipelines.pca_train_pipeline,
#     display_name="Train PCA Model",
#     help_text="Train PCA Model",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="pca-resize-embed",
#     pipeline_func=pipelines.pca_resize_and_embed_pipeline,
#     display_name="Resize PCA model and embed data",
#     help_text="Resize PCA model and embed data",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="pca-resize-full",
#     pipeline_func=pipelines.pca_resize_full_cycle_pipeline,
#     display_name="Resize PCA model and create all resources",
#     help_text="Resize PCA model and create all resources",
# )
#
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="benchmarking-generate-outputs",
#     pipeline_func=pipelines.generate_cas_outputs_pipeline,
#     display_name="Generate CAS Outputs for Benchmarking",
#     help_text="Generate CAS Outputs for Benchmarking",
# )
#
# kubeflow_command_registry.register_pipeline_as_command(
#     typer_app=typer_app,
#     name="benchmarking-calculate-metrics",
#     pipeline_func=pipelines.calculate_metrics_pipeline,
#     display_name="benchmarking_calculate_metrics",
#     help_text="Benchmark CAS results",
# )
#
#
#
# if __name__ == "__main__":
#     typer_app()


cli = click.Group()

# Register the pipeline command
kubeflow_command_registry.register_pipeline_as_command(
    click_group=cli,
    name="bq-ops-prepare-and-extract",
    pipeline_func=pipelines.bq_ops_prepare_and_extract,
    display_name="BQ Ops prepare and extract",
    help_text="BQ Ops prepare extract tables and process extract",
)

# Run the CLI
if __name__ == "__main__":
    cli()
