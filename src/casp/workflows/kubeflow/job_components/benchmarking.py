from kfp import dsl


@dsl.component(
    packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@0.0.3"],
)
def benchmark_cas(gcs_config_path: str) -> None:
    """
    Run the benchmarking for CAS.

    :param gcs_config_path: Config file path on GCS.
    """
    import pickle

    import yaml
    from smart_open import open

    from casp.scripts.benchmarking import benchmark_cas
    from casp.workflows.kubeflow.job_components_library import clients, matching_engine_client

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    cas_api_token = config_data["cas_api_token"]
    cas_api_url = config_data["cas_api_url"]
    wandb_project = config_data["wandb_project"]
    benchmarking_data_paths = config_data["benchmark_data_paths"].split(",")
    model_name = config_data["model_name"]
    cell_ontology_resource_path = config_data["cell_ontology_resource_path"]
    cl_labels_to_names_map_path = config_data["cl_labels_to_names_map_path"]
    metrics_metadata_path = config_data["metrics_metadata_path"]

    # Arguments for deploying and undeploying index
    index_region = config_data["index_region"]
    index_name = config_data["index_name"]
    project_id = config_data["project_id"]
    index_endpoint_id = config_data["index_endpoint_id"]
    index_id = config_data["index_id"]
    display_name = config_data["index_display_name"]
    deployed_index_id = config_data["deployed_index_id"]

    matching_engine_client.CustomMatchingEngineIndex.deploy_index(
        region=index_region,
        project_id=project_id,
        index_id=index_id,
        index_endpoint_id=index_endpoint_id,
        display_name=display_name,
        deployed_index_id=deployed_index_id,
    )
    clients.CellariumCloudInternalServiceClient.update_index_deployed_id(
        index_name=index_name, deployed_index_id=deployed_index_id
    )

    benchmark_cas(
        benchmarking_data_paths=benchmarking_data_paths,
        model_name=model_name,
        cas_api_token=cas_api_token,
        cas_api_url=cas_api_url,
        co_resource_path=cell_ontology_resource_path,
        cl_labels_to_names_map_path=cl_labels_to_names_map_path,
        output_path=metrics_metadata_path,
        wandb_project=wandb_project,
    )

    matching_engine_client.CustomMatchingEngineIndex.undeploy_index(
        region=index_region,
        project_id=project_id,
        index_endpoint_id=index_endpoint_id,
        deployed_index_id=deployed_index_id,
    )
    clients.CellariumCloudInternalServiceClient.update_index_deployed_id(index_name=index_name, deployed_index_id=None)
