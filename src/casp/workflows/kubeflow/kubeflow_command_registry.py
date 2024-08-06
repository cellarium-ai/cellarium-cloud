import typing as t
import tempfile

import typer
from typing_extensions import Annotated
import kfp
from kfp import dsl
from io import StringIO
import yaml

from casp.workflows.kubeflow import config_management
from casp.workflows.kubeflow.pipelines.utils import submit_pipeline


import ast


class DSLParser(ast.NodeVisitor):
    def __init__(self):
        self.components = {}

    def visit_FunctionDef(self, node):
        for decorator in node.decorator_list:
            if isinstance(decorator, ast.Attribute) and decorator.attr == "component":
                component_name = node.name
                self.components[component_name] = []
                self.current_component = component_name
                self.generic_visit(node)
                self.current_component = None
            elif (
                isinstance(decorator, ast.Call)
                and isinstance(decorator.func, ast.Attribute)
                and decorator.func.attr == "component"
            ):
                component_name = node.name
                self.components[component_name] = []
                self.current_component = component_name
                self.generic_visit(node)
                self.current_component = None
        self.generic_visit(node)

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and self.current_component:
            self.components[self.current_component].append(node.func.id)
        self.generic_visit(node)


# Function to compile the pipeline to a temporary YAML file and extract components
def extract_pipeline_data(pipeline_func):
    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as temp_file:
        temp_file_path = temp_file.name
        # Compile the pipeline to a temporary YAML file
        kfp.compiler.Compiler().compile(pipeline_func, temp_file_path)

    # Load the YAML file as a dictionary
    with open(temp_file_path, "r") as file:
        pipeline_dict = yaml.safe_load(file)

    pipeline_data = extract_pipelines_and_components(pipeline_dict)

    # Ensure nested pipelines' components are included correctly
    def clean_pipeline_data(data):
        cleaned_data = {}
        for key, value in data.items():
            if isinstance(value, dict):
                cleaned_data[key] = clean_pipeline_data(value)
            else:
                cleaned_data[key] = list(value)
        return cleaned_data

    return clean_pipeline_data(pipeline_data)


def register_pipeline_as_command(
    typer_app: typer.Typer, name: str, pipeline_func: t.Callable, display_name: str, help_text: str
):
    def command(
        config_yaml_path: Annotated[
            t.Optional[str],
            typer.Option(help="Path to the local YAML config file containing a list of configs for each run."),
        ] = None,
        generate_config_template: Annotated[
            bool, typer.Option("--template", help="If applied, config template will be generated")
        ] = False,
    ):
        components = extract_pipeline_data(pipeline_func=pipeline_func)
        print("ASDASDASDAS)", components)
        if generate_config_template:
            print("printing template")

        elif config_yaml_path:
            print("EXECUTING!")
            # config_paths = config_management.create_configs_pipelines_as_components(config_yaml_path)
            # submit_pipeline(
            #     pipeline_func=pipeline_func,
            #     pipeline_display_name=display_name,
            #     pipeline_kwargs={"pipeline_config_paths": config_paths},
            # )

    typer_app.command(name=name, help=help_text)(command)
