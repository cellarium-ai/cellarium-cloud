import argparse

import fastavro
from google.api_core.exceptions import Conflict
from google.cloud import bigquery
from google.cloud import storage
import os


def create_table(client, project, dataset, tablename, schema, clustering_fields):
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
    missing = list(filter(lambda f: not os.path.exists(f), filenames))
    if len(missing) > 0:
        raise ValueError(
            f"Missing Avro files for loading to BigQuery: {', '.join(missing)}")


def bucket_and_prefix(project, gcs_prefix):
    import re
    client = storage.Client(project=project)
    (bucket_name, object_prefix) = re.search(r'gs://([^/]+)/(.*)$', gcs_prefix).groups()
    bucket = client.bucket(bucket_name)
    return bucket, object_prefix.rstrip('/')


def create_bigquery_objects(client, project, dataset):
    create_dataset(client, project, dataset, "US")

    create_table(client, project, dataset, "cas_load_status",
                 [
                     bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("cas_table_name", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("load_timestamp", "TIMESTAMP", mode="NULLABLE"),
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_ingest_info",
                 [
                     bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED"),
                     # TODO get direct JSON metadata loading working
                     # bigquery.SchemaField("uns_metadata", "JSON", mode="REQUIRED"),
                     bigquery.SchemaField("uns_metadata", "STRING", mode="REQUIRED"),
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_cell_info",
                 [
                     bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("original_cell_id", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("cell_type", "STRING", mode="REQUIRED"),
                     # TODO get direct JSON metadata loading working
                     # bigquery.SchemaField("obs_metadata", "JSON", mode="REQUIRED"),
                     bigquery.SchemaField("obs_metadata", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED")
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_feature_info",
                 [
                     bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("original_feature_id", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("feature_name", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("var_metadata", "STRING", mode="REQUIRED"),
                     # TODO get direct JSON metadata loading working
                     # bigquery.SchemaField("var_metadata", "JSON", mode="REQUIRED"),
                     bigquery.SchemaField("cas_ingest_id", "STRING", mode="REQUIRED")
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_raw_count_matrix",
                 [
                     bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("raw_counts", "INTEGER", mode="REQUIRED")
                 ],
                 ["cas_cell_index"]
                 )


def process(project, dataset, avro_prefix, gcs_prefix, force_bq_append):
    client = bigquery.Client(project=project)
    create_bigquery_objects(client, project, dataset)

    input_file_types = ['ingest_info', 'cell_info', 'feature_info', 'raw_counts']
    input_filenames = [f'{avro_prefix}_{file_type}.avro' for file_type in input_file_types]
    confirm_input_files_exist(input_filenames)
    ingest_filename, cell_filename, feature_filename, raw_counts_filename = input_filenames

    # Grab the `cas_ingest_id` from the ingest file.
    with open(ingest_filename, 'rb') as f:
        reader = fastavro.reader(f)
        ingest_id = next(reader)['cas_ingest_id']

    (bucket, object_prefix) = bucket_and_prefix(project, gcs_prefix)

    def stage_file(file_to_stage):
        print(f"Staging '{file}' to '{gcs_prefix}/{file}'...")
        blob = bucket.blob(f'{object_prefix}/{file_to_stage}')
        blob.upload_from_filename(file_to_stage)
        print(f"Staged '{file}'.")

    def unstage_file(file_to_unstage):
        blob = bucket.blob(f'{object_prefix}/{file_to_unstage}')
        blob.delete()
        print(f"Removed staged file '{gcs_prefix}/{file}'.")

    staged_files = []
    gcs_prefix = gcs_prefix.rstrip('/')
    pairs = [
        ("ingest_info", ingest_filename),
        ("cell_info", cell_filename),
        ("feature_info", feature_filename),
        ("raw_count_matrix", raw_counts_filename)
    ]
    job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.AVRO)

    for table, file in pairs:
        table = f"cas_{table}"
        table_id = f'{project}.{dataset}.{table}'
        load_table_id = f'{project}.{dataset}.cas_load_status'

        # Check if this data needs to be loaded
        # noinspection SqlResolve
        print(f"Checking load status of table '{table_id}'...")
        query = f"""SELECT COUNT(*) AS row_count FROM `{load_table_id}` """ + \
                f"""WHERE cas_ingest_id = "{ingest_id}" AND cas_table_name = "{table}" """
        job = client.query(query)
        row_count = 0
        for row in job.result():
            row_count = row['row_count']

        if row_count > 0:
            print(f"'{table_id}' was already loaded, skipping.")
            continue

        stage_file(file)
        staged_files.append(file)

        uri = f'{gcs_prefix}/{file}'
        load_job = client.load_table_from_uri(
            uri, table_id, job_config=job_config
        )  # Make an API request.
        result = load_job.result()  # Waits for the job to complete.
        print(f"{result.output_rows} rows loaded into BigQuery table '{table_id}'.")

        # Record the load of this data
        # noinspection SqlResolve
        print(f"Recording '{table_id}' as loaded")
        query = \
            f"""INSERT INTO `{project}.{dataset}.cas_load_status`(cas_ingest_id, cas_table_name, load_timestamp)""" + \
            f""" VALUES ("{ingest_id}", "{table}", CURRENT_TIMESTAMP()) """
        job = client.query(query)
        job.result()

    for f in staged_files:
        unstage_file(f)

    print("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Initialize CASP tables')

    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--avro_prefix', type=str, help='Prefix with which Avro files are named', required=True)
    parser.add_argument('--gcs_prefix', type=str, help='GCS prefix to which Avro files should be staged', required=True)
    parser.add_argument('--force_bq_append', type=bool,
                        help='Append data to BigQuery tables even if data some data is already loaded', required=False)

    args = parser.parse_args()

    process(args.project, args.dataset, args.avro_prefix, args.gcs_prefix, args.force_bq_append)
