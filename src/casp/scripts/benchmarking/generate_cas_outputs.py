import pickle
import typing as t

import anndata
from cellarium.cas import CASClient
from smart_open import open

from casp.scripts.benchmarking import utils


def generate_cas_outputs(
    dataset_paths: t.List[str],
    model_name: str,
    cas_api_token: str,
    cas_api_url: str,
    cas_results_output_path: str,
):
    """
    Generate CAS responses for given datasets using the specified model and CAS API.
    Results will be saved as pickle files in `{cas_results_output_path}/{model_name}/` directory.

    :param dataset_paths: List of paths to the dataset files.
    :param model_name: Name of the model to use for CAS.
    :param cas_api_token: API token for authentication with the CAS API.
    :param cas_api_url: URL of the CAS API.
    :param cas_results_output_path: Path to save the CAS results.
    """
    _dataset_file_paths = utils.get_paths(paths=dataset_paths)
    for dataset_file_path in _dataset_file_paths:
        print(f"Running CAS over {dataset_file_path} dataset")

        cas_client = CASClient(api_token=cas_api_token, api_url=cas_api_url)

        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
            matrix=adata,
            cas_model_name=model_name,
            chunk_size=1000,
        )
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        output_path = f"{cas_results_output_path}/{model_name}/cas_output_{dataset_file_name}.pickle"

        with open(output_path, "wb") as f:
            pickle.dump(cas_result, f)
