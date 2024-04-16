import typer
from typing_extensions import Annotated

from casp.bq_scripts import extract_bins_in_parallel_workers, precalculate_fields, prepare_extract

typer_app = typer.Typer()


@typer_app.command()
def bq_ops_precalculate_fields(
    project_id: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    fields: Annotated[str, typer.Option()],
) -> None:
    """
    Precalculate fields in BigQuery.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param fields: Comma separated list of fields in a single string.
    """
    fields = fields.split(",")
    precalculate_fields(project=project_id, dataset=dataset, fields=fields)


@typer_app.command()
def bq_ops_prepare_extract(
    project_id: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    extract_table_prefix: Annotated[str, typer.Option()],
    fq_allowed_original_feature_ids: Annotated[str, typer.Option()],
    extract_bin_size: Annotated[int, typer.Option()],
    filters_json_path: Annotated[str, typer.Option()],
    obs_columns_to_include: Annotated[str, typer.Option()],
    bucket_name: Annotated[str, typer.Option()],
    extract_bucket_path: Annotated[str, typer.Option()],
    ci_random_seed_offset: Annotated[int, typer.Option()] = 0,
    ci_partition_bin_count: Annotated[int, typer.Option()] = 40000,
    ci_partition_size: Annotated[int, typer.Option()] = 10,
) -> None:
    """
    Prepare CAS BigQuery tables used for data extraction.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param extract_table_prefix: A prefix string for naming tables or for similar purposes.
    :param fq_allowed_original_feature_ids: BigQuery table with feature schema needed for the extract
    :param extract_bin_size: The size for the bins where cells are allocated.
    :param bucket_name: GCS Bucket name where to store the metadata files.
    :param extract_bucket_path: GCS Bucket path where the extract will be executed. Used to save metadata files there.
        It is required to use the same bucket path during extract
    :param ci_random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization
        Used in cas_cell_info table.
        `Defaults:` ``0``
    :param ci_partition_bin_count: The count of bins for partitioning cas_cell_info table.
        `Default:` ``40000``
    :param ci_partition_size: The size for the partitions in cas_cell_info table.
        `Default:` ``10``
    :param filters_json_path: A path to json with filters. Filters that have to be included in a SQL query.
        :func:`casp.bq_manager.sql.mako_helpers.parse_where_body`
        Filter format should follow convention described in
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    obs_columns_to_include = obs_columns_to_include.split(",")
    prepare_extract(
        project_id=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        filters_json_path=filters_json_path,
        obs_columns_to_include=obs_columns_to_include,
        bucket_name=bucket_name,
        extract_bucket_path=extract_bucket_path,
        ci_random_seed_offset=ci_random_seed_offset,
        ci_partition_bin_count=ci_partition_bin_count,
        ci_partition_size=ci_partition_size,
    )


@typer_app.command()
def bq_ops_extract_bins_in_parallel_workers(
    project_id: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    extract_table_prefix: Annotated[str, typer.Option()],
    start_bin: Annotated[int, typer.Option()],
    end_bin: Annotated[int, typer.Option()],
    output_bucket_name: Annotated[str, typer.Option()],
    extract_bucket_path: Annotated[str, typer.Option()],
    obs_columns_to_include: Annotated[str, typer.Option()],
) -> None:
    """
    Extract bins in parallel workers with using python :class:`concurrent.futures.ThreadPoolExecutor`.

    :param project_id: Google Cloud Project id
    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param start_bin: Start bin to extract
    :param end_bin: End bin to extract
    :param output_bucket_name: Name of GCS bucket
    :param extract_bucket_path: Path where the extract files and subdirectories should be located.
        Should correspond to the directory provided to prepare_extract script as current script uses `shared_meta`
        files.
    :param obs_columns_to_include:  Comma separated columns from `cas_cell_info` table to include in ``adata.obs``.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
            E.g. ``"c"`` for `cas_cell_info`, ``"i"`` for `cas_ingest_info.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    obs_columns_to_include = obs_columns_to_include.split(",")
    extract_bins_in_parallel_workers(
        project_id=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=start_bin,
        end_bin=end_bin,
        output_bucket_name=output_bucket_name,
        extract_bucket_path=extract_bucket_path,
        obs_columns_to_include=obs_columns_to_include,
    )


if __name__ == "__main__":
    typer_app()
