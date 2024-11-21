import click
import typing as t
import importlib
from casp.workflows.kubeflow import config_management
from casp.workflows.kubeflow.pipelines.utils import submit_pipeline


def lazy_import_function(module_name, function_name):
    module = importlib.import_module(module_name)
    return getattr(module, function_name)


def register_pipeline_as_command(
    click_group: click.Group, name: str, pipeline_func_name: str, display_name: str, help_text: str
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
        pipeline_func = lazy_import_function(
            module_name="casp.workflows.kubeflow.pipelines.pipelines", function_name=pipeline_func_name
        )
        submit_pipeline(
            pipeline_func=pipeline_func,
            pipeline_display_name=display_name,
            pipeline_kwargs=config_paths,
        )
