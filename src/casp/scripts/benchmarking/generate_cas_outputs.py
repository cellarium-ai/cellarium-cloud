import logging
import pickle
import typing as t

import anndata
from smart_open import open

from casp.scripts.benchmarking import utils

logger = logging.getLogger(__name__)


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
    # Import it here, because this dependency is only needed in this function and only should be installed in the
    # environment where this function is called.
    from cellarium.cas import CASClient

    _dataset_file_paths = utils.get_paths(paths=dataset_paths)
    gcs_output_path = cas_results_output_path.replace("gs://", "")
    gcs_bucket_name = gcs_output_path.split("/")[0]
    gcs_output_path = "/".join(gcs_output_path.split("/")[1:])
    existing_output_files = set(
        f"gs://{gcs_bucket_name}/{x}"
        for x in utils.list_files_in_bucket(bucket_name=gcs_bucket_name, prefix=gcs_output_path)
    )
    for dataset_file_path in _dataset_file_paths:
        logger.info(f"Running CAS over {dataset_file_path} dataset...")
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        output_path = f"{cas_results_output_path}/{model_name}/cas_output_{dataset_file_name}.pickle"
        if output_path in existing_output_files:
            logger.info(f"Skipping dataset {dataset_file_name} as its output already exists")
            continue
        cas_client = CASClient(api_token=cas_api_token, api_url=cas_api_url)

        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        cas_result = cas_client.annotate_matrix_cell_type_ontology_aware_strategy(
            matrix=adata,
            cas_model_name=model_name,
            chunk_size=1000,
        )

        logger.info("Uploading output to bucket...")

        with open(output_path, "wb") as f:
            pickle.dump(cas_result, f)
