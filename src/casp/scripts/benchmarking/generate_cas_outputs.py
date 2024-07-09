import pickle
import typing as t
import time

import anndata
from cellarium.cas import CASClient, exceptions as cas_exceptions
from smart_open import open

from casp.scripts.benchmarking import utils


class AnnotationError(Exception):
    pass


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
    delay = 30
    num_attempts = 5

    for dataset_file_path in _dataset_file_paths:
        print(f"Running CAS over {dataset_file_path} dataset...")

        cas_client = CASClient(api_token=cas_api_token, api_url=cas_api_url)

        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        curr_attempt = 0
        done = False
        while curr_attempt <= num_attempts and not done:
            try:
                cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
                    matrix=adata,
                    cas_model_name=model_name,
                    chunk_size=1000,
                )
            except cas_exceptions.DataProcessingError:
                time.sleep(delay)
                curr_attempt += 1
                if curr_attempt < num_attempts:
                    print(f"Retrying annotation, number of current retry: {curr_attempt}...")
                    continue

                raise AnnotationError("Exceeded number of allowed retries. Exiting.")
            else:
                done = True
                print("Successfully finished annotation.")

        print("Uploading output to bucket...")
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        output_path = f"{cas_results_output_path}/{model_name}/cas_output_{dataset_file_name}.pickle"

        with open(output_path, "wb") as f:
            pickle.dump(cas_result, f)
