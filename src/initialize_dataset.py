from google.api_core.exceptions import Conflict
from google.cloud import bigquery
import argparse


def create_table(client, project, dataset, tablename, schema, clustering_fields):
    table_id = f"{project}.{dataset}.{tablename}"
    
    table = bigquery.Table(table_id, schema=schema)
    if clustering_fields:
        table.clustering_fields = clustering_fields

    try:
        table = client.create_table(table)  # Make an API request.
        print(f"Created {table_id}")
    except Conflict:
        print(f"Table {table_id} exists, continuing")


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


def process(project, dataset):
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Initialize CASP tables')

    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)

    # Execute the parse_args() method
    args = parser.parse_args()

    process(args.project, args.dataset)
