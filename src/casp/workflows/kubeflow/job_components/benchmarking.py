from kfp import dsl

from casp.workflows.kubeflow import machine_specs


@dsl.component(
    base_image=machine_specs.DOCKER_IMAGE_NAME_CPU,
    packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@1.4.1-alpha.3"],
)
def benchmark_cas(gcs_config_path: str) -> None:
    """
    Run the benchmarking for CAS.

    :param gcs_config_path: Config file path on GCS.
    """
    import anndata
    import neptune
    import pandas as pd
    import yaml
    from cellarium.cas import CASClient
    from smart_open import open

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    neptune_project = config_data["neptune_project"]
    neptune_api_token = config_data["neptune_api_token"]
    benchmarking_data_paths = config_data["benchmark_data_paths"].split(",")
    model_name = config_data["model_name"]
    cas_api_token = config_data["cas_api_token"]

    run = neptune.init_run(
        project=neptune_project,
        api_token=neptune_api_token,
        name=f"benchmark-{model_name}",
        tags=["benchmarking", model_name],
    )
    run["parameters"]
    cas_client = CASClient(api_token=cas_api_token)

    result_df = pd.DataFrame(columns=["top_1", "top_3", "top_5"])

    for data_path in benchmarking_data_paths:
        with open(data_path, "rb") as file:
            adata = anndata.read_h5ad(file)

        ground_truths = adata.obs["cell_type"].values
        cas_res = cas_client.annotate_anndata(adata=adata, cas_model_name=model_name, chunk_size=1000)

        for ground_truth, match_result in zip(ground_truths, cas_res):
            matches = match_result["matches"]
            predicted_cell_types = []

            for i, match in enumerate(matches):
                predicted_cell_types.append(match["cell_type"])

            top_1 = 1 if ground_truth in predicted_cell_types[:1] else 0
            top_3 = 1 if ground_truth in predicted_cell_types[:3] else 0
            top_5 = 1 if ground_truth in predicted_cell_types[:5] else 0

            new_row = pd.DataFrame(
                {
                    "top_1": [top_1],
                    "top_3": [top_3],
                    "top_5": [top_5],
                    "ground_truth": [ground_truth],
                    "predicted_cell_types": [", ".join(predicted_cell_types)],
                    "query_cell_id": [match_result["query_cell_id"]],
                }
            )
            result_df = pd.concat([result_df, new_row], ignore_index=True)

    run["benchmarking_result"]["top_1"] = result_df["top_1"].mean()
    run["benchmarking_result"]["top_3"] = result_df["top_3"].mean()
    run["benchmarking_result"]["top_5"] = result_df["top_5"].mean()

    csv_file_name = "result.csv"
    result_df.to_csv(csv_file_name, index=False)
    run["data/benchmarking-result-df"].upload(csv_file_name)
