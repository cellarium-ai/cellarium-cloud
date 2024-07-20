from kfp import dsl


@dsl.component(packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@0.0.7", "owlready2"])
def generate_cas_outputs(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts import benchmarking
    from casp.workflows.kubeflow.job_components_library import clients, matching_engine_client

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    cas_api_token = config_data["cas_api_token"]
    cas_api_url = config_data["cas_api_url"]
    dataset_paths = config_data["dataset_paths"]
    model_name = config_data["model_name"]
    cas_results_output_path = config_data["cas_results_output_path"]

    # Arguments for deploying and undeploying index
    index_region = config_data["index_region"]
    index_name = config_data["index_name"]
    project_id = config_data["project_id"]
    index_endpoint_id = config_data["index_endpoint_id"]
    index_id = config_data["index_id"]
    display_name = config_data["index_display_name"]
    deployed_index_id = config_data["deployed_index_id"]
    index_endpoint_full_id = f"projects/350868384795/locations/us-central1/indexEndpoints/{index_endpoint_id}"

    matching_engine_client.CustomMatchingEngineIndex.deploy_index(
        region=index_region,
        project_id=project_id,
        index_id=index_id,
        index_endpoint_id=index_endpoint_id,
        display_name=display_name,
        deployed_index_id=deployed_index_id,
    )
    clients.CellariumCloudInternalServiceClient.update_index_with_deployment_info(
        index_name=index_name,
        update_kwargs={"deployed_index_id": deployed_index_id, "endpoint_id": index_endpoint_full_id},
    )
    try:
        benchmarking.generate_cas_outputs(
            dataset_paths=dataset_paths,
            model_name=model_name,
            cas_api_token=cas_api_token,
            cas_api_url=cas_api_url,
            cas_results_output_path=cas_results_output_path,
        )
    except Exception as e:
        print(f"Smth went wrong: {e}")
        raise e
    finally:
        matching_engine_client.CustomMatchingEngineIndex.undeploy_index(
            region=index_region,
            project_id=project_id,
            index_endpoint_id=index_endpoint_id,
            deployed_index_id=deployed_index_id,
        )
        clients.CellariumCloudInternalServiceClient.update_index_with_deployment_info(
            index_name=index_name,
            update_kwargs={"deployed_index_id": deployed_index_id},
        )


@dsl.component(packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@0.0.7", "owlready2"])
def calculate_metrics(gcs_config_path: str) -> None:
    """
    Run the benchmarking for CAS.

    :param gcs_config_path: Config file path on GCS.
    """
    import yaml
    from smart_open import open


    from casp.scripts import benchmarking
    from casp.scripts.benchmarking import utils

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    wandb_project = config_data["wandb_project"]
    dataset_paths = config_data["dataset_paths"]
    cas_result_paths = config_data["cas_result_paths"]
    model_name = config_data["model_name"]
    num_hops = config_data["num_hops"]
    cell_ontology_resource_path = config_data["cell_ontology_resource_path"]
    cl_labels_to_names_map_path = config_data["cl_labels_to_names_map_path"]
    metrics_metadata_path = config_data["metrics_metadata_path"]
    batch_size = config_data["batch_size"]

    import wandb
    import pandas as pd

    run = wandb.init(project=wandb_project, name=f"benchmarking_model_{model_name}")

    _dataset_paths = utils.get_paths(paths=dataset_paths)
    _cas_result_paths = utils.get_paths(paths=cas_result_paths)

    metrics_dir = f"{metrics_metadata_path}/{model_name}"
    results = []

    for dataset_file_path, cas_result_path in zip(_dataset_paths, _cas_result_paths):
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        metrics_output_filepath = f"{metrics_dir}/metrics_{dataset_file_name}.csv"
        df = pd.read_csv(metrics_output_filepath)

    # benchmarking.calculate_metrics_for_cas_responses(
    #     dataset_paths=dataset_paths,
    #     cas_result_paths=cas_result_paths,
    #     model_name=model_name,
    #     co_resource_path=cell_ontology_resource_path,
    #     cl_labels_to_names_map_path=cl_labels_to_names_map_path,
    #     output_path=metrics_metadata_path,
    #     wandb_project=wandb_project,
    #     batch_size=batch_size,
    #     num_hops=num_hops,
    # )
