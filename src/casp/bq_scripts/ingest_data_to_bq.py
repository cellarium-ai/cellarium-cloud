import argparse

import fastavro
from google.api_core.exceptions import Conflict
from google.cloud import bigquery, storage


def ingest_data_to_bq(project, dataset, avro_prefix, gcs_stage_dir, credentials=None):
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

    # Grab the `cas_ingest_id` from the ingest file.
    with open(ingest_filename, "rb") as file:
        reader = fastavro.reader(file)
        ingest_id = next(reader)["cas_ingest_id"]

    gcs_stage_dir = gcs_stage_dir.rstrip("/")
    pairs = [
        ("ingest_info", ingest_filename, bigquery.SourceFormat.AVRO),
        ("cell_info", cell_filename, bigquery.SourceFormat.AVRO),
        ("feature_info", feature_filename, bigquery.SourceFormat.AVRO),
        ("raw_count_matrix", raw_counts_pattern, bigquery.SourceFormat.CSV),
    ]

    for table, file_pattern, file_format in pairs:
        table = f"cas_{table}"
        table_id = f"{project}.{dataset}.{table}"

        uri = f"{gcs_stage_dir}/{file_pattern}"
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, description="Initialize CASP tables")

    parser.add_argument("--project", type=str, help="BigQuery Project", required=True)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--avro_prefix", type=str, help="Prefix with which Avro files are named", required=True)
    parser.add_argument(
        "--gcs_path_prefix", type=str, help="GCS prefix to which Avro files should be staged", required=True
    )

    args = parser.parse_args()

    ingest_data_to_bq(args.project, args.dataset, args.avro_prefix, args.gcs_path_prefix)
