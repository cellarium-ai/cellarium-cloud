import argparse
import typing as t

from casp.bq_scripts import prepare_all_cell_types, prepare_extract, prepare_measured_genes_info
from casp.services import utils


def main(
    dataset: str,
    extract_table_prefix: str,
    fq_allowed_original_feature_ids: str,
    extract_bin_size: int,
    bucket_name: str,
    filter_by_organism: t.Optional[str] = None,
    filter_by_datasets: t.Optional[str] = None,
    filter_by_is_primary_data: t.Optional[bool] = None,
    obs_columns_to_include: t.Optional[str] = None,
):
    """
    Prepare extract tables in BigQuery data repository.
    This is just a straight wrapper around `casp.bq_scripts.prepare_extract` with only passing
    credentials from the project settings.

    :param dataset: Working BigQuery Dataset
    :param extract_table_prefix: Prefix for extract tables
    :param bucket_name: Bucket where to store dataset info
    :param fq_allowed_original_feature_ids: Fully qualified reference to table of allowed feature names
    :param extract_bin_size: Desired cells per extract bin
    :param filter_by_organism: Optional filter to specify an organism. If not provided, no filtering occurs.
        `Default:` ``None``
    :param filter_by_datasets: Optional filter to specify datasets by their names. If ``None``, no filtering occurs.
        `Default:` ``None``
    :param filter_by_is_primary_data: Optional filter to determine if data is primary or not.
        `Default:` ``None``
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``. If
        not ``None``, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        `Default:` ``None``
    """
    credentials, project_id = utils.get_google_service_credentials()
    filter_by_datasets_split = None if filter_by_datasets is None else filter_by_datasets.split(",")
    obs_columns_to_include_split = None if obs_columns_to_include is None else obs_columns_to_include.split(",")
    prepare_extract(
        project=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        credentials=credentials,
        filter_by_organism=filter_by_organism,
        filter_by_datasets=filter_by_datasets_split,
        filter_by_is_primary_data=filter_by_is_primary_data,
        obs_columns_to_include=obs_columns_to_include_split,
    )
    measured_genes_info_df = prepare_measured_genes_info(
        project=project_id,
        dataset=dataset,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        credentials=credentials,
    )
    all_cell_types_df = prepare_all_cell_types(project=project_id, dataset=dataset, credentials=credentials)
    measured_genes_file_name = "measured_genes_info.csv"
    measured_genes_info_df.to_csv(measured_genes_file_name)

    all_cell_types_file_name = "all_cell_types.csv"
    all_cell_types_df.to_csv(all_cell_types_file_name, index=False)

    utils.upload_file_to_bucket(
        local_file_name=measured_genes_file_name,
        bucket=bucket_name,
        blob_name=f"{dataset}_{extract_table_prefix}_info/{measured_genes_file_name}",
    )
    utils.upload_file_to_bucket(
        local_file_name=all_cell_types_file_name,
        bucket=bucket_name,
        blob_name=f"{dataset}_{extract_table_prefix}_info/{all_cell_types_file_name}",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Prepare CASP tables ML Training/Inference Extract"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix for extract tables", required=True)
    parser.add_argument("--bucket_name", type=str, help="Bucket Name where to store dataset info files", required=True)
    parser.add_argument(
        "--extract_bin_size", type=int, help="desired cells per extract bin", default=10000, required=False
    )
    parser.add_argument(
        "--fq_allowed_original_feature_ids",
        type=str,
        help="fully qualified reference to table of allowed feature names",
        default="dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A",
        required=False,
    )
    parser.add_argument(
        "--filter_by_organism",
        type=str,
        help="Organism to filter cells by that we want to use for preparing extract tables",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--filter_by_datasets",
        type=str,
        help="Comma separated dataset_filename(s) to use to filter by for preparing extract tables",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--filter_by_is_primary_data",
        type=bool,
        help="Is primary data boolean flag form bigquery tables",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--obs_columns_to_include",
        type=str,
        required=False,
        default=None,
    )

    args = parser.parse_args()
    main(
        dataset=args.dataset,
        extract_table_prefix=args.extract_table_prefix,
        fq_allowed_original_feature_ids=args.fq_allowed_original_feature_ids,
        extract_bin_size=args.extract_bin_size,
        bucket_name=args.bucket_name,
        filter_by_organism=args.filter_by_organism,
        filter_by_datasets=args.filter_by_datasets,
        filter_by_is_primary_data=args.filter_by_is_primary_data,
        obs_columns_to_include=args.obs_columns_to_include,
    )
