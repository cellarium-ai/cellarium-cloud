import argparse
import json
import typing as t

from smart_open import open

from casp.bq_scripts import prepare_all_cell_types, prepare_extract, prepare_measured_genes_info
from casp.services import utils


def main(
    dataset: str,
    extract_table_prefix: str,
    fq_allowed_original_feature_ids: str,
    extract_bin_size: int,
    bucket_name: str,
    extract_bucket_path: str,
    filters_json_path: t.Optional[str] = None,
    obs_columns_to_include: t.Optional[str] = None,
):
    """
    Prepare extract tables in BigQuery data repository.
    This is just a straight wrapper around `casp.bq_scripts.prepare_extract` with only passing
    credentials from the project settings.

    :param dataset: Working BigQuery Dataset
    :param extract_table_prefix: Prefix for extract tables
    :param fq_allowed_original_feature_ids: Fully qualified reference to table of allowed feature names
    :param extract_bin_size: Desired cells per extract bin
    :param bucket_name: Bucket where to store dataset info
    :param extract_bucket_path: Path where the extract files and subdirectories should be located.
        Should correspond to the directory provided to prepare_extract script as current script uses `shared_meta`
        files.
    :param filters_json_path: file path to the JSON containing filter criteria, structured as
        {column_name__filter_type: value}. |br|
        Supported filter_types: |br|
            ``"eq"`` - Used for an 'equals' comparison. |br|
                Example: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``. |br|
            ``"in"`` - Used for an 'in' comparison with a set of values. |br|
                Example: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type in
                    ('T cell', 'neuron')``. |br|
        `Default:` ``None``
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``. If
        not ``None``, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`. |br|
        `Default:` ``None``
    """
    credentials, project_id = utils.get_google_service_credentials()
    obs_columns_to_include_split = None if obs_columns_to_include is None else obs_columns_to_include.split(",")

    filters = None
    if filters_json_path is not None:
        with open(filters_json_path) as f:
            filters = json.loads(f.read())

    prepare_extract(
        project=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        credentials=credentials,
        filters=filters,
        obs_columns_to_include=obs_columns_to_include_split,
    )
    measured_genes_info_df = prepare_measured_genes_info(
        project=project_id,
        dataset=dataset,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        credentials=credentials,
    )
    # all_cell_types_df = prepare_all_cell_types(project=project_id, dataset=dataset, credentials=credentials)
    measured_genes_file_name = "measured_genes_info.csv"
    measured_genes_info_df.to_csv(measured_genes_file_name)

    # all_cell_types_file_name = "all_cell_types.csv"
    # all_cell_types_df.to_csv(all_cell_types_file_name, index=False)

    utils.upload_file_to_bucket(
        local_file_name=measured_genes_file_name,
        bucket=bucket_name,
        blob_name=f"{extract_bucket_path}/shared_meta/{measured_genes_file_name}",
    )
    # utils.upload_file_to_bucket(
    #     local_file_name=all_cell_types_file_name,
    #     bucket=bucket_name,
    #     blob_name=f"{extract_bucket_path}/shared_meta/{all_cell_types_file_name}",
    # )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Prepare CASP tables ML Training/Inference Extract"
    )
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix for extract tables", required=True)
    parser.add_argument("--bucket_name", type=str, help="Bucket Name where to store dataset info files", required=True)
    parser.add_argument("--extract_bucket_path", type=str, required=True)
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
    parser.add_argument("--filters_json_path", type=str, required=False, default=None)
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
        extract_bucket_path=args.extract_bucket_path,
        filters_json_path=args.filters_json_path,
        obs_columns_to_include=args.obs_columns_to_include,
    )
