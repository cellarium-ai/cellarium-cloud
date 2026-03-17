import math
import pathlib
import time

import fastavro
from google.api_core.exceptions import Forbidden
from google.cloud import bigquery
from smart_open import open

from casp.scripts.bq_ops import create_bigquery_objects
from casp.services import utils

MAX_RETRY_ATTEMPTS = 5


def _ingest_data_to_bq(project, dataset, avro_prefix, gcs_bucket_name, gcs_stage_dir, credentials=None):
    """
    Main method that drives the 5 high level steps of BigQuery data loading.
    """

    if credentials is None:
        client = bigquery.Client(project=project)
    else:
        client = bigquery.Client(project=project, credentials=credentials)

    input_file_types = ["ingest_info", "cell_info", "feature_info", "raw_counts"]
    input_filenames = [f"{avro_prefix}_{file_type}.avro" for file_type in input_file_types]
    # TODO: reimplement
    ingest_filename, cell_filename, feature_filename, raw_counts_filename = input_filenames
    raw_counts_pattern = f"{avro_prefix}_raw_counts.*.csv"

    gcs_stage_dir = gcs_stage_dir.rstrip("/")
    # Grab the `cas_ingest_id` from the ingest file.
    utils.download_file_from_bucket(
        bucket_name=gcs_bucket_name,
        source_blob_name=f"{gcs_stage_dir}/{ingest_filename}",
        destination_file_name=ingest_filename,
    )
    with open(ingest_filename, "rb") as file:
        reader = fastavro.reader(file)
        ingest_id = next(reader)["cas_ingest_id"]

    pairs = [
        ("ingest_info", ingest_filename, bigquery.SourceFormat.AVRO),
        ("cell_info", cell_filename, bigquery.SourceFormat.AVRO),
        ("feature_info", feature_filename, bigquery.SourceFormat.AVRO),
        ("raw_count_matrix", raw_counts_pattern, bigquery.SourceFormat.CSV),
    ]

    for table, file_pattern, file_format in pairs:
        table = f"cas_{table}"
        table_id = f"{project}.{dataset}.{table}"

        uri = f"gs://{gcs_bucket_name}/{gcs_stage_dir}/{file_pattern}"
        job_config = bigquery.LoadJobConfig(source_format=file_format)

        load_job = client.load_table_from_uri(uri, table_id, job_config=job_config)  # Make an API request.
        result = load_job.result()  # Waits for the job to complete.
        print(f"{result.output_rows} rows loaded into BigQuery table '{table_id}'.")

    print("Updating ingest timestamp")
    # noinspection SqlResolve
    query = f"""

    UPDATE `{dataset}.cas_ingest_info` SET ingest_timestamp = CURRENT_TIMESTAMP()
        WHERE cas_ingest_id = "{ingest_id}"

    """
    job = client.query(query)
    job.result()

    print("Done.")


def ingest_data_to_bq(
    project_id: str,
    gcs_bucket_name: str,
    dataset: str,
    gcs_stage_dir: str,
    max_retry_attempts: int = MAX_RETRY_ATTEMPTS,
):
    """
    Ingest files prepared by `bq_ops.anndata_to_avro` script. If error happens during ingest, the script would retry
    ``max_retry_attempts`` times.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param gcs_bucket_name: GCS Bucket name
    :param dataset: BigQuery dataset name
    :param gcs_stage_dir: GCS directory where ingest files are stored
    :param max_retry_attempts: Maximum number for retries per ingest |br|
        `Default:` ``5``
    """
    bq_client = bigquery.Client()
    create_bigquery_objects(client=bq_client, project=project_id, dataset=dataset)

    ingest_file_blobs = utils.list_blobs(bucket_name=gcs_bucket_name, prefix=gcs_stage_dir)
    blob_names = [x.name for x in ingest_file_blobs]
    ingest_avro_prefixes = []

    for blob_name in blob_names:
        path = pathlib.Path(blob_name)
        blob_directory = str(path.parent)

        # Getting prefixes (everything before `_cell_info`, `_ingest_info`, `_feature_info` `_raw_counts`)
        suffixes = ["_cell_info", "_ingest_info", "_feature_info", "_raw_counts"]
        prefix = path.name

        for suffix in suffixes:
            prefix = prefix.split(suffix)[0]

        if blob_directory != gcs_stage_dir:
            continue

        ingest_avro_prefixes.append(prefix)

    ingest_avro_prefixes = list(set(ingest_avro_prefixes))

    for avro_prefix in ingest_avro_prefixes:
        need_retry = True
        attempt_counter = 1

        while need_retry and attempt_counter <= max_retry_attempts:
            try:
                print(f"Ingesting files with prefix: {avro_prefix}")
                _ingest_data_to_bq(
                    project=project_id,
                    dataset=dataset,
                    gcs_bucket_name=gcs_bucket_name,
                    avro_prefix=avro_prefix,
                    gcs_stage_dir=gcs_stage_dir,
                )
            except Forbidden as e:
                # It can happen when limit of number of table update operations is exceeded
                # Retrying the operation 5 times with an exponential backoff as suggested in the docs:
                # https://cloud.google.com/bigquery/quotas#standard_tables
                print("Was not able to ingest data", e)
                time_to_wait = math.exp(attempt_counter)
                time.sleep(time_to_wait)
                if attempt_counter <= 5:
                    print("Retrying another attempt...")
            except Exception:
                need_retry = False
            else:
                need_retry = False

            attempt_counter += 1
