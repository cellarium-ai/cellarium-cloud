from kfp import dsl
import yaml
from casp.workflows.kubeflow import machine_specs


@dsl.component(
    base_image=machine_specs.DOCKER_IMAGE_NAME_CPU,
    packages_to_install=["git+https://github.com/cellarium-ai/cellarium-cas.git@1.4.1rc1"]
)
def benchmark_cas(gcs_config_path: str) -> None:
    """
    Run the benchmarking for CAS.

    :param gcs_config_path: Config file path on GCS.
    """
    import neptune
    from smart_open import open
    import anndata
    import pandas as pd
    from cellarium.cas import CASClient

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    neptune_project = config_data["neptune_project"]
    neptune_api_token = config_data["neptune_api_token"]
    benchmarking_data_paths = config_data["benchmarking_data_paths"]
    model_name = config_data["model_name"]
    cas_api_token = config_data["cas_api_token"]

    run = neptune.init_run(
        project=neptune_project,
        api_token=neptune_api_token,
        name="benchmarking",
        tags=["benchmarking", model_name],
    )

    cas_client = CASClient(api_token=cas_api_token)

    for data_path in benchmarking_data_paths:
        with open(data_path, "r") as file:
            adata = anndata.read_h5ad(file)

        ground_truths = adata.obs["cell_type"].values
        cas_res = cas_client.annotate(adata=adata, cas_model_name=model_name, chunk_size=1000)
        result_df = pd.DataFrame(columns=["top_1", "top_3", "top_5"])

        for ground_truth, match_result in zip(ground_truths, cas_res):
            matches = match_result["matches"]
            predicted_cell_types = []

            for i, match in enumerate(matches):
                predicted_cell_types.append(match["cell_type"])

            top_1 = 1 if ground_truth in predicted_cell_types[:1] else 0
            top_3 = 1 if ground_truth in predicted_cell_types[:3] else 0
            top_5 = 1 if ground_truth in predicted_cell_types[:5] else 0

            result_df = result_df.append(
                {
                    "top_1": top_1,
                    "top_3": top_3,
                    "top_5": top_5,
                    "ground_truth": ground_truth,
                    "predicted_cell_types": ", ".join(predicted_cell_types),
                    "query_cell_id": match_result["query_cell_id"],
                },
                ignore_index=True,
            )


