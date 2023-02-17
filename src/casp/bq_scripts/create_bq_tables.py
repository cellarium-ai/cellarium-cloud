"""
* Creates BigQuery dataset if necessary
* Creates tables in BiqQuery dataset if necessary
"""

import os
import re

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
