import pickle
import typing as t

import anndata
import pandas as pd
import wandb
from cellarium.cas import CASClient
from smart_open import open

from casp.scripts.benchmarking import utils


def main(
    benchmarking_data_paths: t.List[str],
    model_name: str,
    cas_api_token: str,
    cas_api_url: str,
    co_resource_path: str,
    cl_labels_to_names_map_path: str,
    output_path: str,
    wandb_project: str,
):
    with open(co_resource_path, "rb") as f:
        co_resource = pickle.load(f)

    with open(cl_labels_to_names_map_path, "rb") as f:
        cl_labels_to_names_map = pickle.load(f)

    cas_client = CASClient(api_token=cas_api_token, api_url=cas_api_url)

    run = wandb.init(project=wandb_project, name=f"benchmarking_model_{model_name}")

    metrics_dir = f"{output_path}/{model_name}"
    results = []
    for dataset_file_path in benchmarking_data_paths:
        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        if adata.raw is not None:
            # Getting raw expression counts instead of normalized.
            adata = adata.raw.to_adata()

        ground_truths = adata.obs.cell_type.values
        cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
            matrix=adata,
            cas_model_name=model_name,
            chunk_size=1000,
        )
        df_metrics = utils.calculate_metrics_for_cas_output(
            ground_truths=ground_truths,
            cas_result=cas_result,
            co_resource=co_resource,
            cl_labels_to_names_map=cl_labels_to_names_map,
            num_hops=4,
        )
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        df_metrics.to_csv(f"{metrics_dir}/metrics_{dataset_file_name}.csv")
        results.append(df_metrics)

    means_list = [df.mean() for df in results]
    means_df = pd.DataFrame(means_list)
    aggregated_means = means_df.mean()
    run.log(aggregated_means.to_dict())
