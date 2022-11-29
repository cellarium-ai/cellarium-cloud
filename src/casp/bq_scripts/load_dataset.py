"""
* Creates BigQuery dataset if necessary
* Creates tables in BiqQuery dataset if necessary
* Stages input Avro files to GCS
* Loads BigQuery tables from Avro files
* Unstages input Avro files
"""

import argparse
import glob
import os
import re

import fastavro
from google.api_core.exceptions import Conflict
from google.cloud import bigquery, storage


def create_table(client, project, dataset, tablename, schema, clustering_fields):
    """
    Create the specified table in the specified project / dataset using the specified schema and clustering fields.
    """
    table_id = f"{project}.{dataset}.{tablename}"

    table = bigquery.Table(table_id, schema=schema)
    if clustering_fields:
        table.clustering_fields = clustering_fields

    try:
        _ = client.create_table(table)  # Make an API request.
        print(f"Created '{table_id}'.")
    except Conflict:
        print(f"Table '{table_id}' exists, continuing.")


def create_dataset(client, project, dataset, location):
    """
    Create the specified dataset in the specified project and location.
    """
    dataset_id = f"{project}.{dataset}"

    # Construct a full Dataset object to send to the API.
    dataset = bigquery.Dataset(dataset_id)
    dataset.location = location

    # Send the dataset to the API for creation, with an explicit timeout.
    # Raises google.api_core.exceptions.Conflict if the Dataset already
    # exists within the project.
    try:
        _ = client.create_dataset(dataset, timeout=30)  # Make an API request.
        print(f"Created dataset {dataset_id}.")
    except Conflict:
        print(f"Dataset {dataset_id} exists, continuing.")


def confirm_input_files_exist(filenames):
    """
    Validate that the specified input files exist or throw.
    """
    missing = list(filter(lambda f: not os.path.exists(f), filenames))
    if len(missing) > 0:
        raise ValueError(f"Missing Avro files for loading to BigQuery: {', '.join(missing)}")


def bucket_and_prefix(project, gcs_path_prefix):
    """
    Extract a GCS bucket and object prefix from the specified GCS bucket + prefix.
    """
    client = storage.Client(project=project)
    match = re.search(r"gs://([^/]+)/(.*)$", gcs_path_prefix)
    if not match:
        raise ValueError(f"Invalid gcs_path_prefix `{gcs_path_prefix}`; should look like `gs://bucket/prefix1/etc`")
    (bucket_name, object_prefix) = match.groups()
    bucket = client.bucket(bucket_name)
    return bucket, object_prefix.rstrip("/")


def create_bigquery_objects(client, project, dataset):
    """
    Create the core BigQuery dataset and tables for CASP.
    """
    create_dataset(client, project, dataset, "US")

    create_table(
        client,
        project,
        dataset,
        "cas_ingest_info",
        [
            bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("uns_metadata", "JSON", mode="REQUIRED"),
            bigquery.SchemaField("ingest_timestamp", "TIMESTAMP", mode="NULLABLE"),
        ],
        [],
    )

    create_table(
        client,
        project,
        dataset,
        "cas_cell_info",
        [
            bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("original_cell_id", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("cell_type", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("obs_metadata", "JSON", mode="REQUIRED"),
            bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED"),
        ],
        [],
    )

    create_table(
        client,
        project,
        dataset,
        "cas_feature_info",
        [
            bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("original_feature_id", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("feature_name", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("var_metadata", "JSON", mode="REQUIRED"),
            bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED"),
        ],
        [],
    )

    create_table(
        client,
        project,
        dataset,
        "cas_raw_count_matrix",
        [
            bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("raw_counts", "INTEGER", mode="REQUIRED"),
        ],
        ["cas_cell_index"],
    )


def process(project, dataset, avro_prefix, gcs_path_prefix):
    """
    Main method that drives the 5 high level steps of BigQuery data loading.
    """
    (bucket, object_prefix) = bucket_and_prefix(project, gcs_path_prefix)

    client = bigquery.Client(project=project)
    create_bigquery_objects(client, project, dataset)

    input_file_types = ["ingest_info", "cell_info", "feature_info", "raw_counts"]
    input_filenames = [f"{avro_prefix}_{file_type}.avro" for file_type in input_file_types]
##    confirm_input_files_exist(input_filenames)
    ingest_filename, cell_filename, feature_filename, raw_counts_filename = input_filenames
    raw_counts_pattern = f"{avro_prefix}_raw_counts.*.csv"

    # Grab the `cas_ingest_id` from the ingest file.
    with open(ingest_filename, "rb") as file:
        reader = fastavro.reader(file)
        ingest_id = next(reader)["cas_ingest_id"]

    def stage_file(file_to_stage):
        print(f"Staging '{file_to_stage}' to '{gcs_path_prefix}/{file_to_stage}'...")
        blob = bucket.blob(f"{object_prefix}/{file_to_stage}")
        blob.upload_from_filename(file_to_stage)
        print(f"Staged '{file_to_stage}'.")

    def unstage_file(file_to_unstage):
        blob = bucket.blob(f"{object_prefix}/{file_to_unstage}")
        blob.delete()
        print(f"Removed staged file '{gcs_path_prefix}/{file_to_unstage}'.")

    staged_files = []
    gcs_path_prefix = gcs_path_prefix.rstrip("/")
    pairs = [
       ("ingest_info", ingest_filename, bigquery.SourceFormat.AVRO),
       ("cell_info", cell_filename, bigquery.SourceFormat.AVRO),
       ("feature_info", feature_filename, bigquery.SourceFormat.AVRO),
        ("raw_count_matrix", raw_counts_pattern, bigquery.SourceFormat.CSV),
    ]
    job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.AVRO)

    for table, file_pattern, format in pairs:
        table = f"cas_{table}"
        table_id = f"{project}.{dataset}.{table}"

        for file in glob.glob(file_pattern):
            stage_file(file)
            staged_files.append(file)

        uri = f"{gcs_path_prefix}/{file_pattern}"
        job_config = bigquery.LoadJobConfig(source_format=format)

        load_job = client.load_table_from_uri(uri, table_id, job_config=job_config)  # Make an API request.
        result = load_job.result()  # Waits for the job to complete.
        print(f"{result.output_rows} rows loaded into BigQuery table '{table_id}'.")

    for file in staged_files:
        unstage_file(file)

    print("Updating ingest timestamp")
    # noinspection SqlResolve
    query = f"""

    UPDATE `{dataset}.cas_ingest_info` SET ingest_timestamp = CURRENT_TIMESTAMP()
        WHERE cas_ingest_id = "{ingest_id}"

    """
    job = client.query(query)
    job.result()

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, description="Initialize CASP tables")

    parser.add_argument("--project", type=str, help="BigQuery Project", required=True)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--avro_prefix", type=str, help="Prefix with which Avro files are named", required=True)
    parser.add_argument(
        "--gcs_path_prefix", type=str, help="GCS prefix to which Avro files should be staged", required=True
    )

    args = parser.parse_args()

    process(args.project, args.dataset, args.avro_prefix, args.gcs_path_prefix)