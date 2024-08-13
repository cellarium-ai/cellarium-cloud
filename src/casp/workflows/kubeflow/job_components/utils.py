import typing as t

from kfp import dsl
from kfp.components import BaseComponent
from kfp.dsl.python_component import PythonComponent
from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component

from casp.workflows.kubeflow import machine_specs_utils, exceptions


def get_dsl_component_from_python_function(
    component_function: t.Callable[..., t.Any], base_image: str
) -> t.Union[t.Callable, PythonComponent, BaseComponent]:
    return dsl.component(base_image=base_image)(component_function)


def create_job(
    component_func: t.Callable[..., t.Any], component_name: str, gcs_config_path: str
) -> t.Callable[[], t.Any]:
    """
    Create a custom training Google Vertex AI job for running a custom training component.

    :param component_func: Custom training component.
    :param component_name: Name of the component. Used to get machine specs.
    :param gcs_config_path: Config file path on GCS.

    :return: Callable custom training job.
    """
    machine_specs_info = machine_specs_utils.read_machine_specs()
    # Use Dummy Specs when `machine_specs_info` is an empty dictionary. This is needed because Kubeflow
    # indexes and executes the pipelines while the pipeline got initialized
    try:
        machine_spec = machine_specs_info[component_name]
    except KeyError:
        print(machine_specs_info)
        raise exceptions.MachineSpecsDoesntExist(f"No machine specs found for component {component_name}")

    dsl_component = get_dsl_component_from_python_function(
        component_function=component_func, base_image=machine_spec["base_image"]
    )

    job = create_custom_training_job_from_component(
        dsl_component,
        display_name=machine_spec["display_name"],
        replica_count=machine_spec["replica_count"],
        machine_type=machine_spec["machine_type"],
        accelerator_type=machine_spec["accelerator_type"],
        accelerator_count=machine_spec["accelerator_count"],
        service_account=machine_spec.get("service_account", ""),
        boot_disk_size_gb=machine_spec.get("boot_disk_size_gb", 100),
        env=machine_spec.get("env", None),
    )

    return lambda: job(gcs_config_path=gcs_config_path)
