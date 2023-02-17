import time
import math
import argparse
from google.cloud import bigquery
from google.api_core.exceptions import Forbidden
from casp.bq_scripts import create_bigquery_objects, ingest_data_to_bq
from casp.services import utils


def get_avro_prefixes(bucket_name, gcs_stage_dir):
    ingest_file_blobs = utils.list_blobs(bucket_name=bucket_name, prefix=gcs_stage_dir)
    blob_names = [x.name.replace(gcs_stage_dir, "") for x in ingest_file_blobs]
    return set(x.split("_")[0].lstrip("/") for x in blob_names)


def main(dataset: str, gcs_bucket_name: str, gcs_stage_dir: str, delete_ingest_files: bool = False):
    credentials, project_id = utils.get_google_service_credentials()
    bq_client = bigquery.Client(project=project_id, credentials=credentials)
    create_bigquery_objects(client=bq_client, project=project_id, dataset=dataset)
    ingest_avro_prefixes = get_avro_prefixes(bucket_name=gcs_bucket_name, gcs_stage_dir=gcs_stage_dir)
    print(ingest_avro_prefixes)
    for avro_prefix in ingest_avro_prefixes:
        need_retry = True
        attempt_counter = 1

        while need_retry or attempt_counter <= 5:
            try:
                ingest_data_to_bq(
                    project=project_id,
                    dataset=dataset,
                    gcs_bucket_name=gcs_bucket_name,
                    avro_prefix=avro_prefix,
                    gcs_stage_dir=gcs_stage_dir,
                    credentials=credentials
                )
            except Forbidden:
                # It can happen when limit of number of table update operations is exceeded
                # Retrying the operation 5 times with an exponential backoff as suggested in the docs:
                # https://cloud.google.com/bigquery/quotas#standard_tables
                print("Was not able to ingest data")
                time_to_wait = math.exp(attempt_counter)
                time.sleep(time_to_wait)
                if attempt_counter <= 5:
                    print("Retrying another attempt...")
            else:
                need_retry = False

            attempt_counter += 1

        if delete_ingest_files:
            print("Deleting ingest files...")
            utils.delete_folder_from_bucket(bucket_name=gcs_bucket_name, folder_name=gcs_stage_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Ingest files to BigQuery from a staged bucket having avro and csv files"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--gcs_bucket_name", type=str, help="A name of working GCS bucket", required=True)
    parser.add_argument(
        "--gcs_stage_dir", type=str, help="A folder name in GCS bucket where all of the ingest files are stored",
        required=True
    )
    parser.add_argument(
        "--delete_ingest_files",
        type=bool,
        help="Whether or not delete ingest files after",
        default=False,
        required=False
    )
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        gcs_bucket_name=args.gcs_bucket_name,
        gcs_stage_dir=args.gcs_stage_dir,
        delete_ingest_files=args.delete_ingest_files
    )
