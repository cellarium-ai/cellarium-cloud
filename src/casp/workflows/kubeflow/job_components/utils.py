import typing as t

from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component

from casp.workflows.kubeflow import machine_specs

# TODO: Think of passing a bare python function in `create_job` and the
# decorating it inside of `create_job` with `@dsl.component`


def create_job(
    dsl_component: t.Callable[..., t.Any], component_name: str, gcs_config_path: str
) -> t.Callable[[], t.Any]:
    """
    Create a custom training Google Vertex AI job for running a custom training component.

    :param dsl_component: Custom training component.
    :param component_name: Name of the component. Used to get machine specs.
    :param gcs_config_path: Config file path on GCS.

    :return: Callable custom training job.
    """
    machine_spec = machine_specs.component_machine_specs_map[component_name]
    # TODO:
    # dsl_component.base_image = machine_spec.base_image
    job = create_custom_training_job_from_component(
        dsl_component,
        display_name=machine_spec.display_name,
        replica_count=machine_spec.replica_count,
        machine_type=machine_spec.machine_type,
        accelerator_type=machine_spec.accelerator_type,
        accelerator_count=machine_spec.accelerator_count,
    )
    return lambda: job(gcs_config_path=gcs_config_path)
