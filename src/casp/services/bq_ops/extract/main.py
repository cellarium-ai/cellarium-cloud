import argparse
import concurrent.futures as concurrency
import logging
import multiprocessing
import os

from casp import bq_scripts
from casp.services import utils


def extract_task(
    dataset: str,
    extract_table_prefix: str,
    bin_number: int,
    file_name: str,
    output_bucket_name: str,
    output_bucket_directory: str,
) -> None:
    """
    Wrapper task `casp.bq_scripts.extract_minibatch_to_anndata` which processes exactly
    one bin at a time and saves the output anndata file to a GCS bucket.

    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param bin_number: Bin to extract
    :param file_name: Name for a local file (.h5ad anndata) to save the output
    :param output_bucket_name: Name of GCS bucket
    :param output_bucket_directory: Directory in the bucket to save outputs
    """
    credentials, project_id = utils.get_google_service_credentials()
    bq_scripts.extract_minibatch_to_anndata(
        project=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=bin_number,
        end_bin=bin_number,
        output=file_name,
        credentials=credentials,
    )
    blob_name = f"{output_bucket_directory}/{file_name}"
    utils.upload_file_to_bucket(local_file_name=file_name, blob_name=blob_name, bucket=output_bucket_name)
    os.remove(file_name)
    logging.info(msg=f"Processed bin {bin_number}")


def main(
    dataset: str,
    extract_table_prefix: str,
    start_bin: int,
    end_bin: int,
    output_bucket_name: str,
    output_bucket_directory: str,
) -> None:
    """
    Extract anndatafiles from bigquery extract tables. Run extract tasks concurrently: 1 task per CPU core
    This is a wrapper around `casp.bq_scripts.extract_minibatch_to_anndata`
    Features:
    1. Concurrently activate `casp.bq_scripts.extract_minibatch_to_anndata` multiple times (wrapped script)
    2. Upload script output to the bucket (wrapper after)

    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param start_bin: Starting (inclusive) integer bin to extract
    :param end_bin: Ending (inclusive) integer bin to extract
    :param output_bucket_name: Name of GCS bucket
    :param output_bucket_directory: Directory in the bucket to save outputs
    """
    num_of_workers = multiprocessing.cpu_count()
    with concurrency.ProcessPoolExecutor(max_workers=num_of_workers) as executor:
        futures = []
        for bin_number in range(start_bin, end_bin + 1, 1):
            file_name = f"extract_{bin_number}.h5ad"
            extract_task_kwargs = {
                "dataset": dataset,
                "extract_table_prefix": extract_table_prefix,
                "bin_number": bin_number,
                "file_name": file_name,
                "output_bucket_name": output_bucket_name,
                "output_bucket_directory": output_bucket_directory,
            }

            future = executor.submit(extract_task, **extract_task_kwargs)
            futures.append(future)

        done, not_done = concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)

        if not_done:
            print("The following futures have not completed: \n", not_done)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, description="Extract Service Execution")
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix of extract tables", required=True)
    parser.add_argument("--start_bin", type=int, help="Starting (inclusive) integer bin to extract", required=True)
    parser.add_argument("--end_bin", type=int, help="Ending (inclusive) integer bin to extract", required=True)
    parser.add_argument("--output_bucket_name", type=str, help="Bucket where to save script results", required=True)
    parser.add_argument(
        "--output_bucket_directory", type=str, help="A specific location in a bucket for the results", required=True
    )
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        extract_table_prefix=args.extract_table_prefix,
        start_bin=args.start_bin,
        end_bin=args.end_bin,
        output_bucket_name=args.output_bucket_name,
        output_bucket_directory=args.output_bucket_directory,
    )
