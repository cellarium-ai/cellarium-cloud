import argparse
import os
import secrets
from casp.bq_scripts import anndata_to_avro
from casp.services import utils


def main(
        dataset: str,
        gcs_bucket_name: str,
        gcs_file_path: str,
        cas_cell_index_start: int,
        cas_feature_index_start: int,
        gcs_ingest_path: str,
) -> None:
    filename = f"{secrets.token_hex()}.h5ad"
    utils.download_file_from_bucket(
        bucket_name=gcs_bucket_name,
        source_blob_name=gcs_file_path,
        destination_file_name=filename
    )
    credentials, project_id = utils.get_google_service_credentials()

    prefix = secrets.token_hex()

    anndata_to_avro(
        input_file=filename,
        project=project_id,
        dataset=dataset,
        cas_cell_index_start=cas_cell_index_start,
        cas_feature_index_start=cas_feature_index_start,
        prefix=prefix,
        load_uns_data=False
    )
    ingest_files = [x for x in os.listdir(os.curdir) if x.startswith(prefix)]

    for filename in ingest_files:
        utils.upload_file_to_bucket(
            local_file_name=filename,
            bucket=gcs_bucket_name,
            blob_name=f"{gcs_ingest_path}/{filename}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Convert AnnData Single Cell Expression Data into format for loading into BQ"
    )

    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=False)
    parser.add_argument("--gcs_input_bucket", type=str, help="AnnData format input file", required=True)
    parser.add_argument("--gcs_file_path", type=str, help="AnnData format input file", required=True)
    parser.add_argument(
        "--gcs_stage_prefix", type=str, help="GCS prefix to which Avro files should be staged", required=True
    )

    parser.add_argument("--cas_cell_index_start", type=int, help="starting number for cell index", required=False)
    parser.add_argument("--cas_feature_index_start", type=int, help="starting number for feature index", required=False)

    args = parser.parse_args()

    main(
        dataset=args.dataset,
        gcs_bucket_name=args.gcs_input_bucket,
        gcs_file_path=args.gcs_file_path,
        cas_cell_index_start=args.cas_cell_index_start,
        cas_feature_index_start=args.cas_feature_index_start,
        gcs_ingest_path=args.gcs_stage_prefix
    )
