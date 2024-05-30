"""Create extract file(s) in a google bucket. To be run locally."""

from google.cloud import bigquery
import argparse
import json
import os
from casp.services import utils

credentials, project = utils.get_google_service_credentials()


def create_table_with_all_ensembl_ids(dataset: str) -> str:
    client = bigquery.Client(project=project, credentials=credentials)

    # create a new dataset called `{dataset}_schemas` if it doesn't exist
    dataset_id = f"{dataset}_schemas"
    dataset_ref = client.dataset(dataset_id)
    if not client.get_dataset(dataset_ref):
        try:
            dataset = client.create_dataset(dataset_ref)
            print(f"Dataset {dataset.dataset_id} created.")
        except Exception as e:
            print(f"Error creating dataset `{dataset_id}`: {str(e)}")

    # create a table with all ensembl ids
    table_name = f"{dataset}_schemas.gene_schema"
    query = f"""CREATE OR REPLACE TABLE `{table_name}` AS
        SELECT
            original_feature_id AS feature_name,
            ROW_NUMBER() OVER () AS index
        FROM (
            SELECT DISTINCT original_feature_id
            FROM `{dataset}.cas_feature_info`
        ) AS distinct_features
        ORDER BY index;
    """
    client.query(query)
    return table_name


def get_max_bin_number(dataset: str, extract_table_prefix: str):
    client = bigquery.Client(project=project, credentials=credentials)
    query = f"SELECT MAX(extract_bin) FROM `{dataset}.{extract_table_prefix}__extract_cell_info`"
    result = client.query(query).result()
    df = result.to_dataframe()
    max_extract_bin = df.iloc[0].item()
    return max_extract_bin


def main(
    dataset: str,
    extract_table_prefix: str,
    extract_bucket_name: str,
    extract_bucket_path: str, 
    extract_bin_size: int,
    filters_json_path: str,
    obs_columns_to_include: str,
    max_processes: int,
    cleanup: bool,
):

    this_filepath = os.path.abspath(__file__)
    prepare_extract_script = os.path.join(os.path.dirname(this_filepath), "..", "..", "services/bq_ops/prepare_extract/main.py")
    extract_script = os.path.join(os.path.dirname(this_filepath), "..", "..", "services/bq_ops/extract/main.py")
    # fq_allowed_original_feature_ids = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A"
    try:
        fq_allowed_original_feature_ids = project + "." + create_table_with_all_ensembl_ids(dataset=dataset)
    except Exception as e:
        print("ERROR!", str(e))
        print(f"It is likely that you need to manually create the dataset `{dataset}_schemas` using the BigQuery UI")
        return

    # run preparation
    prepare_extract_cmd = f"""python {prepare_extract_script} 
        --dataset {dataset} 
        --extract_table_prefix {extract_table_prefix} 
        --bucket_name {extract_bucket_name} 
        --extract_bucket_path {extract_bucket_path}
        --extract_bin_size {extract_bin_size}
        --fq_allowed_original_feature_ids {fq_allowed_original_feature_ids} 
        --filters_json_path {filters_json_path} 
        --obs_columns_to_include {obs_columns_to_include}
    """
    print(prepare_extract_cmd)
    os.system(prepare_extract_cmd.replace("\n", " "))

    print('Preparation complete. Beginning h5ad creation.')

    # run extraction
    start_bin = 0
    end_bin = get_max_bin_number(dataset=dataset, extract_table_prefix=extract_table_prefix)
    extract_cmd = f"""python {extract_script} 
        --dataset {dataset} 
        --extract_table_prefix {extract_table_prefix} 
        --start_bin {start_bin} 
        --end_bin {end_bin} 
        --output_bucket_name {extract_bucket_name} 
        --extract_bucket_path {extract_bucket_path} 
        --obs_columns_to_include {obs_columns_to_include}
        --max_processes {max_processes}
    """
    print(f"Extracting bins {start_bin} to {end_bin}")
    print(extract_cmd)
    os.system(extract_cmd.replace("\n", " "))

    if cleanup:
        print('Cleaning up intermediate tables...')
        client = bigquery.Client(project=project, credentials=credentials)
        for suffix in ["__extract_cell_info", "__extract_feature_info", "__extract_raw_count_matrix"]:
            table_name = f"{dataset}.{extract_table_prefix}{suffix}"
            client.delete_table(table_name, not_found_ok=True)
            print(f'Intermediate table `{table_name}` deleted.')

    print('Finished data extraction.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Extract (filtered) data from BigQuery to a bucket as h5ad file(s)"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix for new BigQuery Tables related to this extract", required=True)
    parser.add_argument("--gcs_bucket_name", type=str, help="Base GCS bucket name for extracted h5ad files, no gs://", required=True)
    parser.add_argument(
        "--gcs_h5ad_dir",
        type=str,
        help="A folder name in the GCS bucket where all of the output h5ad files will be stored",
        required=True,
    )
    parser.add_argument("--filters_json", type=str, help="JSON file defining the obs filters", required=True)
    parser.add_argument("--obs_columns_to_include", type=str, help="Comma-separated list (no spaces) of obs column names", required=True)
    parser.add_argument("--extract_bin_size", type=int, default=10_000, help="Number of cells in each h5ad file")
    parser.add_argument("--max_processes", type=int, default=4, help="Maximum number of concurrent processes to run")
    parser.add_argument("--cleanup", action='store_true', help="Include flag to delete intermediate tables after extraction")
    
    args = parser.parse_args()

    # validation
    if "cas_cell_index" in args.obs_columns_to_include:
        raise ValueError("cas_cell_index is always included in the extract, so do not include it in --obs_columns_to_include")
    # check that there is a period in each json field in filters_json
    with open(args.filters_json) as f:
        filters_json = json.load(f)
    for f in filters_json:
        if "." not in f:
            raise ValueError(f"Filters json ({args.filters_json}) field '{f}' must be fully qualified: if this is an obs column, use 'c.{f}'")

    main(
        dataset=args.dataset,
        extract_table_prefix=args.extract_table_prefix,
        extract_bucket_name=args.gcs_bucket_name,
        extract_bucket_path=args.gcs_h5ad_dir, 
        extract_bin_size=args.extract_bin_size,
        filters_json_path=args.filters_json,
        obs_columns_to_include=args.obs_columns_to_include,
        max_processes=args.max_processes,
        cleanup=args.cleanup,
    )
