import typing as t

import typer
from typing_extensions import Annotated
from casp.workflows.kubeflow import config_management
from casp.workflows.kubeflow.pipelines.utils import submit_pipeline


def register_pipeline_as_command(
        typer_app: typer.Typer, name: str, pipeline_func: t.Callable, display_name: str, help_text: str
):
    def command(
            config_yaml_path: Annotated[
                t.Optional[str],
                typer.Option(help="Path to the local YAML config file containing a list of configs for each run."),
            ] = None,
            # generate_config_template: Annotated[
            #     bool, typer.Option("--template", help="If applied, config template will be generated")
            # ] = False,
    ):
        config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
        submit_pipeline(
            pipeline_func=pipeline_func,
            pipeline_display_name=display_name,
            pipeline_kwargs=config_paths,
        )

    typer_app.command(name=name, help=help_text)(command)
