import argparse
from google.api_core.exceptions import Conflict
from google.cloud import bigquery
import subprocess


def create_table(client, project, dataset, tablename, schema, clustering_fields):
    table_id = f"{project}.{dataset}.{tablename}"
    
    table = bigquery.Table(table_id, schema=schema)
    if clustering_fields:
        table.clustering_fields = clustering_fields

    try:
        table = client.create_table(table)  # Make an API request.
        print(f"Created '{table_id}'")
    except Conflict:
        print(f"Table '{table_id}' exists, continuing")


def create_dataset(client, project, dataset, location):
    dataset_id = f"{project}.{dataset}"

    # Construct a full Dataset object to send to the API.
    dataset = bigquery.Dataset(dataset_id)
    dataset.location = location

    # Send the dataset to the API for creation, with an explicit timeout.
    # Raises google.api_core.exceptions.Conflict if the Dataset already
    # exists within the project.
    try:
        dataset = client.create_dataset(dataset, timeout=30)  # Make an API request.
        print(f"Created dataset {project}.{dataset}")
    except Conflict:
        print(f"Dataset {project}.{dataset} exists, continuing")


def check_avro_files(avro_prefix):
    file_types = ['cell_info', 'feature_info', 'raw_counts']
    filenames = [f'{avro_prefix}_{file_type}.avro' for file_type in file_types]
    import os
    missing = list(filter(lambda f: not os.path.exists(f), filenames))
    if len(missing) > 0:
        raise ValueError(
            f"Missing Avro files for loading to BigQuery: {', '.join(missing)}")
    return filenames


def process(project, dataset, avro_prefix, gcs_prefix, force_bq_append):
    # Construct a BigQuery client object.
    client = bigquery.Client(project=project)
    
    create_dataset(client, project, dataset, "US")
    
    create_table(client, project, dataset, "cas_cell_info",
        [
            bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("original_cell_id", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("cell_type", "STRING", mode="REQUIRED")
        ],
        []
    )

    create_table(client, project, dataset, "cas_feature_info",
        [
            bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
            bigquery.SchemaField("original_feature_id", "STRING", mode="REQUIRED"),
            bigquery.SchemaField("feature_name", "STRING", mode="REQUIRED")
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

    (cell_filename, feature_filename, raw_counts_filename) = check_avro_files(avro_prefix)
    query = f"""select table_name, total_rows from `{dataset}.INFORMATION_SCHEMA.PARTITIONS` where total_rows > 0"""
    job = client.query(query)
    tables_with_data = [r[0] for r in list(job.result())]

    def stage_file(file_to_stage):
        proc = subprocess.Popen(['gsutil', '-m', 'cp', file_to_stage, gcs_prefix],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        staged_files.append(file_to_stage)

    staged_files = []
    bqload_template = ["bq", "load", "-project_id", project, "--source_format=AVRO"]
    gcs_prefix = gcs_prefix.rstrip('/')
    pairs = [("cell_info", cell_filename), ("feature_info", feature_filename), ("raw_count_matrix", raw_counts_filename)]
    for table, file in pairs:
        table = f"cas_{table}"

        if table in tables_with_data and not force_bq_append:
            table_id = f'{project}.{dataset}.{table}'
            print(f"Table '{table_id}' contains data and `--force_bq_append` not specified, skipping load of '{file}'")
            continue

        stage_file(file)
        print(f"Staged '{file}' to '{gcs_prefix}'")
        proc = subprocess.Popen(
            bqload_template + [f'{dataset}.{table}', f'{gcs_prefix}/{file}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print(f"Loaded '{file}' to '{table}' table")

    if len(staged_files) > 0:
        files = [f'{gcs_prefix}/{f}' for f in staged_files]
        proc = subprocess.Popen(['gsutil', 'rm'] + files,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        print("Removed Avro files from GCS staging area")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Initialize CASP tables')

    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--avro_prefix', type=str, help='Prefix with which Avro files are named', required=True)
    parser.add_argument('--gcs_prefix', type=str, help='GCS prefix to which Avro files should be staged', required=True)
    parser.add_argument('--force_bq_append', type=str,
                        help='Append data to BigQuery tables even if data some data is already loaded', required=False)

    # Execute the parse_args() method
    args = parser.parse_args()

    process(args.project, args.dataset, args.avro_prefix, args.gcs_prefix, args.force_bq_append)
