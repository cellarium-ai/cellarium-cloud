import argparse
import os
import typing as t

from casp.bq_scripts import anndata_to_avro
from casp.services import utils


def main(
    gcs_bucket_name: str,
    gcs_file_path: str,
    cas_cell_index_start: int,
    cas_feature_index_start: int,
    gcs_stage_dir: str,
    original_feature_id_lookup: str,
    load_uns_data: bool = False,
    uns_meta_keys: t.Optional[str] = None,
) -> None:
    """
    Create ingest files that can be read by BigQuery for data ingestion.
    This is a wrapper around `casp.bq_scripts.anndata_to_avro`
    Features:
    1. Download anndata file from the specified bucket location (wrapper before logic)
    2. Activate `casp.bq_scripts.anndata_to_avro` script (wrapped script)
    3. Upload all outputs from the script to the specified bucket (wrapper after logic)

    :param gcs_bucket_name: Working bucket name
    :param gcs_file_path: File path that needs to be processed
    :param cas_cell_index_start:  Starting number for cell index
    :param cas_feature_index_start:  Starting number for feature index
    :param gcs_stage_dir: Name of the directory in google bucket
    :param original_feature_id_lookup: A column name in var dataframe from where to get original feature ids.
        In most of the cases it will be a column with ENSEMBL gene IDs. Default is `index` which means that
        an index column of var dataframe would be used.
    :param load_uns_data: Whether to dump uns
    """
    filename = gcs_file_path.split("/")[-1]
    utils.download_file_from_bucket(
        bucket_name=gcs_bucket_name, source_blob_name=gcs_file_path, destination_file_name=filename
    )
    credentials, project_id = utils.get_google_service_credentials()
    if uns_meta_keys is not None:
        uns_meta_keys_list = set([x.strip() for x in uns_meta_keys.split(",")])
    else:
        uns_meta_keys_list = None

    prefix = filename.split(".")[0]

    anndata_to_avro(
        input_file=filename,
        project=project_id,
        cas_cell_index_start=cas_cell_index_start,
        cas_feature_index_start=cas_feature_index_start,
        prefix=prefix,
        dataset=None,
        load_uns_data=load_uns_data,
        original_feature_id_lookup=original_feature_id_lookup,
        included_adata_uns_keys=uns_meta_keys_list,
    )
    ingest_files = [x for x in os.listdir(os.curdir) if x.startswith(prefix) and not x.endswith(".h5ad")]

    for filename in ingest_files:
        utils.upload_file_to_bucket(
            local_file_name=filename, bucket=gcs_bucket_name, blob_name=f"{gcs_stage_dir}/{filename}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Convert AnnData Single Cell Expression Data into format for loading into BQ"
    )
    parser.add_argument("--gcs_input_bucket", type=str, help="AnnData format input file", required=True)
    parser.add_argument("--gcs_file_path", type=str, help="AnnData format input file", required=True)
    parser.add_argument(
        "--gcs_stage_dir", type=str, help="GCS prefix to which Avro files should be staged", required=True
    )

    parser.add_argument("--cas_cell_index_start", type=int, help="starting number for cell index", required=False)
    parser.add_argument("--cas_feature_index_start", type=int, help="starting number for feature index", required=False)
    parser.add_argument(
        "--original_feature_id_lookup", type=str, help="Column name in `var` for original feature id", required=False
    )
    parser.add_argument(
        "--load_uns_data", type=bool, help="Whether or not to dump uns data", required=False, default=False
    )
    parser.add_argument(
        "--uns_meta_keys", type=str, help="Comma separated set of keys to dump in uns metadata", default=None
    )

    args = parser.parse_args()

    main(
        gcs_bucket_name=args.gcs_input_bucket,
        gcs_file_path=args.gcs_file_path,
        cas_cell_index_start=args.cas_cell_index_start,
        cas_feature_index_start=args.cas_feature_index_start,
        gcs_stage_dir=args.gcs_stage_dir,
        original_feature_id_lookup=args.original_feature_id_lookup,
        load_uns_data=args.load_uns_data,
        uns_meta_keys=args.uns_meta_keys,
    )
