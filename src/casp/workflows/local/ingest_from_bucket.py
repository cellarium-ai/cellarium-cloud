import argparse

from casp.services import utils
import os
import tempfile
from google.cloud import bigquery


def get_max_value(table_name: str, column_name: str) -> float:
    client = bigquery.Client()
    query = f"SELECT MAX({column_name}) FROM {table_name}"
    result = client.query(query).result()
    max_value = result[0][0]
    return max_value


def get_bucket_h5ad_blob_names(bucket_name: str, gcs_ingest_dir: str):
    ingest_h5ad_blobs = utils.list_blobs(bucket_name=bucket_name, prefix=gcs_ingest_dir)
    blob_names = [x.name for x in ingest_h5ad_blobs]
    return blob_names


def main(bucket_name: str, h5ad_file_dir: str, avro_dir: str, dataset: str):

    # blob_names = get_bucket_h5ad_blob_names(bucket_name, h5ad_file_dir)
    with tempfile.NamedTemporaryFile() as f:
        os.system(f"gsutil ls gs://{bucket_name}/{h5ad_file_dir} > {f.name}")
        blob_names = [x.strip() for x in open(f.name).readlines()]
        blob_names = [x for x in blob_names if x.endswith(".h5ad")]

    for blob_name in blob_names:

        print(f"Processing {blob_name}")

        this_filepath = os.path.abspath(__file__)
        avro_script = os.path.join(os.path.dirname(this_filepath), "..", "..", "services/bq_ops/anndata_to_ingest_files/main.py")
        ingest_script = os.path.join(os.path.dirname(this_filepath), "..", "..", "services/bq_ops/ingest_files_to_bq/main.py")
        h5ad_file = blob_name.lstrip(f"gs://{bucket_name}/")

        # run avro creation
        avro_cmd = f"""python {avro_script} 
            --gcs_input_bucket {bucket_name} 
            --gcs_file_path {h5ad_file} 
            --gcs_stage_dir {avro_dir} 
            --dataset {dataset}
            --original_feature_id_lookup index 
            --load_uns_data true 
            --uns_meta_keys 'batch_condition,title,schema_reference,schema_version'
        """
        print(avro_cmd)
        os.system(avro_cmd.replace("\n", " "))

        # run bq ingest
        ingest_cmd = f"""python {ingest_script} 
            --dataset {dataset} 
            --gcs_bucket_name {bucket_name} 
            --gcs_stage_dir {avro_dir} 
            --delete_ingest_files true
        """
        print(ingest_cmd)
        os.system(ingest_cmd.replace("\n", " "))

        print('Finished all files.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Ingest to BigQuery from a bucket of h5ad files following the same schema"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--gcs_bucket_name", type=str, help="Base GCS bucket name, no gs://", required=True)
    parser.add_argument(
        "--gcs_h5ad_dir",
        type=str,
        help="A folder name in GCS bucket where all of the input h5ad files are stored",
        required=True,
    )
    parser.add_argument(
        "--gcs_avro_dir",
        type=str,
        help="A folder name in GCS bucket where all of the temporary avro files are stored",
        required=True,
    )
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        bucket_name=args.gcs_bucket_name,
        h5ad_file_dir=args.gcs_h5ad_dir,
        avro_dir=args.gcs_avro_dir,
    )
