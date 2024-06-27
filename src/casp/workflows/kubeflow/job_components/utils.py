import typing as t

from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component

from casp.workflows.kubeflow import machine_specs_utils

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
    machine_specs_info = machine_specs_utils.read_machine_specs()
    # Use Dummy Specs when `machine_specs_info` is an empty dictionary. This is needed because Kubeflow
    # indexes and executes the pipelines while the pipeline got initialized
    machine_spec = machine_specs_info.get(component_name, machine_specs_utils.DUMMY_SPEC)

    dsl_component.component_spec.implementation.container.image = machine_spec["base_image"]

    # if "pca_index_create" in component_name:
    #     service_account = "fgrab@broadinstitute.org"
    # else:
    service_account = "backend-service-account@dsp-cell-annotation-service.iam.gserviceaccount.com"

    job = create_custom_training_job_from_component(
        dsl_component,
        display_name=machine_spec["display_name"],
        replica_count=machine_spec["replica_count"],
        machine_type=machine_spec["machine_type"],
        accelerator_type=machine_spec["accelerator_type"],
        accelerator_count=machine_spec["accelerator_count"],
        # service_account=machine_spec.get("service_account", ""),
        service_account=service_account,
        boot_disk_size_gb=machine_spec.get("boot_disk_size_gb", 100),
        env=machine_spec.get("env", None)
    )
    return lambda: job(gcs_config_path=gcs_config_path)
