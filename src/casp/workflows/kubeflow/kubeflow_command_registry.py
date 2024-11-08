import click
import typing as t
from casp.workflows.kubeflow import config_management
from casp.workflows.kubeflow.pipelines.utils import submit_pipeline


def register_pipeline_as_command(
    click_group: click.Group, name: str, pipeline_func: t.Callable, display_name: str, help_text: str
):
    @click_group.command(name=name, help=help_text)
    @click.option(
        "--config-yaml-path",
        type=str,
        default=None,
        help="Path to the local YAML config file containing a list of configs for each run.",
    )
    def command(config_yaml_path):
        config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
        submit_pipeline(
            pipeline_func=pipeline_func,
            pipeline_display_name=display_name,
            pipeline_kwargs=config_paths,
        )
