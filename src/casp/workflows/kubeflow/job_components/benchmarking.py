from kfp import dsl


# @dsl.component(
#     packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@0.0.3"],
# )
# def benchmark_cas(gcs_config_path: str) -> None:
#     """
#     Run the benchmarking for CAS.
#
#     :param gcs_config_path: Config file path on GCS.
#     """
#     import anndata
#     import wandb
#     import pandas as pd
#     import yaml
#     from cellarium.cas import CASClient
#     from smart_open import open
#
#     with open(gcs_config_path, "r") as file:
#         config_data = yaml.safe_load(file)
#
#     wandb_project = config_data["wandb_project"]
#     benchmarking_data_paths = config_data["benchmark_data_paths"].split(",")
#     model_name = config_data["model_name"]
#     cas_api_token = config_data["cas_api_token"]
#
#     cas_client = CASClient(api_token=cas_api_token)
#
#     result_df = pd.DataFrame(columns=["top_1", "top_3", "top_5"])
#
#     for data_path in benchmarking_data_paths:
#         with open(data_path, "rb") as file:
#             adata = anndata.read_h5ad(file)
#
#         ground_truths = adata.obs["cell_type"].values
#         cas_res = cas_client.annotate_anndata(adata=adata, cas_model_name=model_name, chunk_size=1000)
#
#         for ground_truth, match_result in zip(ground_truths, cas_res):
#             matches = match_result["matches"]
#             predicted_cell_types = []
#
#             for i, match in enumerate(matches):
#                 predicted_cell_types.append(match["cell_type"])
#
#             top_1 = 1 if ground_truth in predicted_cell_types[:1] else 0
#             top_3 = 1 if ground_truth in predicted_cell_types[:3] else 0
#             top_5 = 1 if ground_truth in predicted_cell_types[:5] else 0
#
#             new_row = pd.DataFrame(
#                 {
#                     "top_1": [top_1],
#                     "top_3": [top_3],
#                     "top_5": [top_5],
#                     "ground_truth": [ground_truth],
#                     "predicted_cell_types": [", ".join(predicted_cell_types)],
#                     "query_cell_id": [match_result["query_cell_id"]],
#                 }
#             )
#             result_df = pd.concat([result_df, new_row], ignore_index=True)
#
#     run["benchmarking_result"]["top_1"] = result_df["top_1"].mean()
#     run["benchmarking_result"]["top_3"] = result_df["top_3"].mean()
#     run["benchmarking_result"]["top_5"] = result_df["top_5"].mean()
#
#     csv_file_name = "result.csv"
#     result_df.to_csv(csv_file_name, index=False)
#     run["data/benchmarking-result-df"].upload(csv_file_name)


@dsl.component(
    packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@0.0.3"],
)
def benchmark_cas(gcs_config_path: str) -> None:
    """
    Run the benchmarking for CAS.

    :param gcs_config_path: Config file path on GCS.
    """
    import pickle
    import anndata
    import wandb
    import pandas as pd
    import yaml
    from cellarium.cas import CASClient
    from smart_open import open
    from google.cloud import aiplatform
    from google.cloud.aiplatform_v1 import IndexEndpointServiceClient
    from google.cloud.aiplatform_v1.types import DeployedIndex, DeployIndexRequest, UndeployIndexRequest

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    wandb_project = config_data["wandb_project"]
    benchmarking_data_paths = config_data["benchmark_data_paths"].split(",")
    model_name = config_data["model_name"]
    cas_api_token = config_data["cas_api_token"]
    cas_api_url = config_data["cas_api_url"]
    index_region = config_data["index_region"]
    project_id = config_data["project_id"]
    index_id = config_data["index_id"]
    index_endpoint_id = config_data["index_endpoint_id"]
    display_name = config_data["index_display_name"]
    deployed_index_id = config_data["deployed_index_id"]

    CELL_ONTOLOGY_RESOURCE_PATH = "benchmarking_resource.pickle"
    CL_LABELS_TO_NAMES_MAP_PATH = "cl_labels_to_names_map.pickle"

    with open(CELL_ONTOLOGY_RESOURCE_PATH, "rb") as f:
        co_resource = pickle.load(f)

    with open(CL_LABELS_TO_NAMES_MAP_PATH, "rb") as f:
        cl_labels_to_names_map = pickle.load(f)

    def deploy_index(region, project_id, index_id, index_endpoint_id, display_name, deployed_index_id):
        # Initialize Vertex AI SDK
        aiplatform.init(
            project=project_id,
            location=region,
        )

        # Define variables
        index_endpoint_name = f"projects/{project_id}/locations/{region}/indexEndpoints/{index_endpoint_id}"
        index_name = f"projects/{project_id}/locations/{region}/indexes/{index_id}"

        # Create client
        client = IndexEndpointServiceClient(client_options={"api_endpoint": f"{region}-aiplatform.googleapis.com"})

        # Create deployed index
        deployed_index = DeployedIndex(
            id=deployed_index_id,
            index=index_name,
            display_name=display_name,
            automatic_resources={"min_replica_count": 1, "max_replica_count": 4},
        )

        # Prepare deployment request
        deploy_index_request = DeployIndexRequest(
            index_endpoint=index_endpoint_name,
            deployed_index=deployed_index,
        )

        print("Starting index deployment...")
        operation = client.deploy_index(request=deploy_index_request)

        # Wait for the operation to complete
        print("Waiting for the deployment to complete...")
        operation.result(timeout=7200)  # Timeout after 2 hours

        # Check if the operation was successful
        if operation.done():
            print("Index deployment completed successfully.")
        else:
            print("Index deployment failed. Please check the Google Cloud Console for details.")

    def benchmark_dataset(adata, model_name):
        cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
            matrix=adata,
            cas_model_name=model_name,
        )

        ground_truths = adata.obs.cell_type.values
        df_result = pd.DataFrame()

        for query_res_obj, ground_truth in zip(cas_result, ground_truths):
            hop_0_sens = hop_1_sens = hop_2_sens = hop_3_sens = 0
            hop_0_spec = hop_1_spec = hop_2_spec = hop_3_spec = 0

            ground_truth_cl_name = cl_labels_to_names_map[ground_truth]

            query_cell_id = query_res_obj["query_cell_id"]

            match_co_data = co_resource[ground_truth_cl_name]
            hop_0 = match_co_data["top_0_ancestors"]
            hop_1 = match_co_data["top_1_ancestors"]
            hop_2 = match_co_data["top_2_ancestors"]
            hop_3 = match_co_data["top_3_ancestors"]

            for match in query_res_obj["matches"]:
                match_cl_name = match["cell_type_ontology_term_id"]
                match_score = match["score"]

                # Hop 0
                if match_cl_name in hop_0:
                    hop_0_sens = max(match_score, hop_0_sens)
                else:
                    hop_0_spec = max(1 - match_score, hop_0_spec)

                # Hop 1
                if match_cl_name in hop_1:
                    hop_1_sens = max(match_score, hop_1_sens)
                else:
                    hop_1_spec = max(match_score, hop_1_spec)

                # Hop 2
                if match_cl_name in hop_2:
                    hop_2_sens = max(match_score, hop_2_sens)
                else:
                    hop_2_spec = max(match_score, hop_2_spec)

                # Hop 3
                if match_cl_name in hop_3:
                    hop_3_sens = max(match_score, hop_3_sens)
                else:
                    hop_3_spec = max(match_score, hop_3_spec)

            hop_0_spec = 1 - hop_0_spec
            hop_1_spec = 1 - hop_1_spec
            hop_2_spec = 1 - hop_2_spec
            hop_3_spec = 1 - hop_3_spec

            df_result = pd.concat(
                [
                    df_result,
                    pd.DataFrame(
                        [
                            {
                                "query_cell_id": query_cell_id,
                                "hop_0_sensitivity": hop_0_sens,
                                "hop_1_sensitivity": hop_1_sens,
                                "hop_2_sensitivity": hop_2_sens,
                                "hop_3_sensitivity": hop_3_sens,
                                "hop_0_specificity": hop_0_spec,
                                "hop_1_specificity": hop_1_spec,
                                "hop_2_specificity": hop_2_spec,
                                "hop_3_specificity": hop_3_spec,
                            }
                        ]
                    ),
                ]
            )

        df_result = df_result.set_index("query_cell_id")
        return df_result

    def undeploy_index(region, project_id, index_endpoint_id, deployed_index_id):
        # Create client
        client = IndexEndpointServiceClient(client_options={"api_endpoint": f"{region}-aiplatform.googleapis.com"})
        index_endpoint_name = f"projects/{project_id}/locations/{region}/indexEndpoints/{index_endpoint_id}"
        # Prepare undeployment request
        undeploy_index_request = UndeployIndexRequest(
            index_endpoint=index_endpoint_name,
            deployed_index_id=deployed_index_id,
        )

        # Undeploy the index and wait for the operation to complete
        print("Starting index undeployment...")
        operation = client.undeploy_index(request=undeploy_index_request)

    deploy_index(
        region=index_region,
        project_id=project_id,
        index_id=index_id,
        index_endpoint_id=index_endpoint_id,
        display_name=display_name,
        deployed_index_id=deployed_index_id,
    )
    cas_client = CASClient(api_token=cas_api_token, cas_api_url=cas_api_url)

    for dataset_file_path in benchmarking_data_paths:
        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        if adata.raw is not None:
            # Getting raw expression counts instead of normalized.
            adata = adata.raw.to_adata()

            df_benchmarking_result = benchmark_dataset(adata=adata, model_name=model_name)

        # cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
        #     matrix=adata,
        #     cas_model_name=model_name,
        # )

    result_df = pd.DataFrame(columns=["top_1", "top_3", "top_5"])
