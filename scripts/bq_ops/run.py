import typer
from typing_extensions import Annotated

from casp.scripts import bq_ops

typer_app = typer.Typer()


@typer_app.command()
def bq_ops_create_ingest_files(
    gcs_bucket_name: Annotated[str, typer.Option()],
    gcs_file_path: Annotated[str, typer.Option()],
    uns_meta_keys: Annotated[str, typer.Option()],
    cas_cell_index_start: Annotated[int, typer.Option()],
    cas_feature_index_start: Annotated[int, typer.Option()],
    load_uns_data: Annotated[bool, typer.Option()],
    original_feature_id_lookup: Annotated[str, typer.Option()],
    dataset_id: Annotated[str, typer.Option()],
    dataset_version_id: Annotated[str, typer.Option()],
    gcs_stage_dir: Annotated[str, typer.Option()],
):
    """
    Create ingest files and upload them in a stage GCS bucket directory. High level entry point, reads the input
    AnnData file and generates Avro files for ingest, cells, features, and raw / core data.

    :param gcs_bucket_name: GCS Bucket name
    :param gcs_file_path: GCS Bucket input file path to process
    :param load_uns_data: Whether to load uns (unstructured) metadata
    :param cas_cell_index_start: Starting number for cell index. If ``None``, increment of maximum index would be used.
        ``None`` requires ``dataset`` to be set as it will use dataset to get the maximum index in the cell_info table
    :param cas_feature_index_start: Starting number for feature index. If ``None``, increment of maximum index would be
        used. ``None`` requires ``dataset`` to be set as it will use dataset to get the maximum index in the cell_info
        table
    :param uns_meta_keys: Comma separated list with a set of keys that need to be dumped in ingest. If None, dump all.
    :param original_feature_id_lookup: A column name in var dataframe from where to get original feature ids.
        In most of the cases it will be a column with ENSEMBL gene IDs. Default is `index` which means that
        an index column of var dataframe would be used.
    :param dataset_id: CZI Dataset ID
    :param dataset_version_id: CZI Dataset version ID
    :param gcs_stage_dir: Stage directory in GCS Bucket where to upload the files

    Example usage
        To create ingest files and upload them to a GCS bucket, use the following command:

        .. code-block:: console

            python casp/scripts/bq_ops/run.py bq-ops-create-ingest-files \\
                --gcs-bucket-name my-bucket \\
                --gcs-file-path path/to/input/file \\
                --uns-meta-keys key1,key2,key3 \\
                --cas-cell-index-start 1000 \\
                --cas-feature-index-start 2000 \\
                --load-uns-data True \\
                --original-feature-id-lookup gene_id \\
                --dataset-id 12345 \\
                --dataset-version-id 67890 \\
                --gcs-stage-dir path/to/stage/dir
    """
    if uns_meta_keys is not None:
        uns_meta_keys = set([x.strip() for x in uns_meta_keys.split(",")])
    else:
        uns_meta_keys = None

    bq_ops.create_ingest_files(
        gcs_bucket_name=gcs_bucket_name,
        gcs_file_path=gcs_file_path,
        uns_meta_keys=uns_meta_keys,
        cas_cell_index_start=cas_cell_index_start,
        cas_feature_index_start=cas_feature_index_start,
        load_uns_data=load_uns_data,
        original_feature_id_lookup=original_feature_id_lookup,
        dataset_id=dataset_id,
        dataset_version_id=dataset_version_id,
        gcs_stage_dir=gcs_stage_dir,
    )


@typer_app.command()
def bq_ops_ingest_data_to_bq(
    project_id: Annotated[str, typer.Option()],
    gcs_bucket_name: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    gcs_stage_dir: Annotated[str, typer.Option()],
    max_retry_attempts: Annotated[int, typer.Option()] = 5,
):
    """
    Ingest files prepared by `bq_ops.anndata_to_avro` script. If error happens during ingest, the script would retry
    ``max_retry_attempts`` times.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param gcs_bucket_name: GCS Bucket name
    :param dataset: BigQuery dataset name
    :param gcs_stage_dir: GCS directory where ingest files are stored
    :param max_retry_attempts: Maximum number for retries per ingest |br|
        `Default:` ``5``

    Example usage
        To ingest data into BigQuery with retries, use the following command:

        .. code-block:: console

            python casp/scripts/bq_ops/run.py bq-ops-ingest-data-to-bq \\
                --project-id my-gcp-project \\
                --gcs-bucket-name my_bucket \\
                --dataset my_dataset \\
                --gcs-stage-dir path/to/ingest/files \\
                --max-retry-attempts 5
    """
    bq_ops.ingest_data_to_bq(
        project_id=project_id,
        gcs_bucket_name=gcs_bucket_name,
        dataset=dataset,
        gcs_stage_dir=gcs_stage_dir,
        max_retry_attempts=max_retry_attempts,
    )


@typer_app.command()
def bq_ops_precalculate_fields(
    project_id: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    fields: Annotated[str, typer.Option()],
) -> None:
    """
    Precalculate fields in BigQuery. You can find more details on which particular fields can be precalculated in
    ``casp/scripts/bq_ops/precalculate_fields.py``.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param fields: Comma separated list of fields in a single string.

    Example usage
        To precalculate a field named "total_mrna_umis" in a BigQuery dataset, use the following command:

        .. code-block:: console

            python casp/scripts/bq_ops/run.py bq-ops-precalculate-fields \\
                --project-id my-gcp-project \\
                --dataset my_dataset \\
                --fields total_mrna_umis
    """
    fields = fields.split(",")
    bq_ops.precalculate_fields(project=project_id, dataset=dataset, fields=fields)


@typer_app.command()
def bq_ops_prepare_extract(
    project_id: Annotated[str, typer.Option()],
    dataset: Annotated[str, typer.Option()],
    extract_table_prefix: Annotated[str, typer.Option()],
    fq_allowed_original_feature_ids: Annotated[str, typer.Option()],
    filters_json_path: Annotated[str, typer.Option()],
    obs_columns_to_include: Annotated[str, typer.Option()],
    bucket_name: Annotated[str, typer.Option()],
    extract_bucket_path: Annotated[str, typer.Option()],
    extract_bin_size: Annotated[int, typer.Option()] = 10000,
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
    :param bucket_name: GCS Bucket name where to store the metadata files.
    :param extract_bucket_path: GCS Bucket path where the extract will be executed. Used to save metadata files there.
        It is required to use the same bucket path during extract
    :param extract_bin_size: The size for the bins where cells are allocated. |br|
        `Default`: ``10000``
    :param ci_random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization
        Used in cas_cell_info table. |br|
        `Default:` ``0``
    :param ci_partition_bin_count: The count of bins for partitioning cas_cell_info table. |br|
        `Default:` ``40000``
    :param ci_partition_size: The size for the partitions in cas_cell_info table. |br|
        `Default:` ``10``
    :param filters_json_path: A path to a JSON file containing filters that should be included in a SQL query.
        The JSON should represent a dictionary containing filter criteria, structured as ``{column_name__filter_type: value}``.
        Supported filter types:

        - ``"eq"``: Used for an 'equals' comparison.

          Example: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``.

        - ``"in"``: Used for an 'in' comparison with a set of values.

          Example: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type IN ('T cell', 'neuron')``.

        - ``"not_eq"``: Used for a 'not equals' comparison, meaning that the query would exclude rows with the specified value.

          Example: ``{"assay__not_eq": "Drop-seq"}`` results in ``assay!='Drop-seq'``.

        - ``"not_in"``: Used for a 'not in' comparison with a set of values to exclude.

          Example: ``{"assay__not_in": ["Drop-seq", "microwell-seq", "BD Rhapsody Targeted mRNA"]}`` results in ``assay NOT IN ('Drop-seq', 'microwell-seq', 'BD Rhapsody Targeted mRNA')``.

        - ``"gt"``: Used for a 'greater than' comparison.

          Example: ``{"total_mrna_umis__gt": 13000}`` results in ``total_mrna_umis > 13000``.

        - ``"gte"``: Used for a 'greater than or equal' comparison.

          Example: ``{"total_mrna_umis__gte": 13000}`` results in ``total_mrna_umis >= 13000``.

        - ``"lt"``: Used for a 'less than' comparison.

          Example: ``{"total_mrna_umis__lt": 13000}`` results in ``total_mrna_umis < 13000``.

        - ``"lte"``: Used for a 'less than or equal' comparison.

          Example: ``{"total_mrna_umis__lte": 13000}`` results in ``total_mrna_umis <= 13000``.

    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``

    Example filters JSON
        The following is an example JSON dictionary for the `filters_json_path` parameter, demonstrating a few filter types:

        .. code-block:: json

            {
                "organism__eq": "Homo sapiens",
                "cell_type__in": ["T cell", "neuron"],
                "assay__not_eq": "Drop-seq",
                "total_mrna_umis__gt": 13000
            }

    Example usage
        To prepare CAS BigQuery tables for data extraction, use the following command:

        .. code-block:: console

            python casp/scripts/bq_ops/run.py bq-ops-prepare-extract \\
                --project-id my-gcp-project \\
                --dataset my_dataset \\
                --extract-table-prefix my_prefix \\
                --fq-allowed-original-feature-ids my_feature_ids_table \\
                --extract-bin-size 10000 \\
                --filters-json-path path/to/filters.json \\
                --obs-columns-to-include "c.cell_type,c.donor_id,c.sex,i.dataset_id" \\
                --bucket-name my_bucket \\
                --extract-bucket-path path/to/extract \\
                --ci-random-seed-offset 0 \\
                --ci-partition-bin-count 40000 \\
                --ci-partition-size 10
    """
    bq_ops.prepare_extract(
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

    Example usage
        To extract bins in parallel workers, use the following command:

        .. code-block:: console

            python casp/scripts/bq_ops/run.py bq-ops-extract-bins-in-parallel-workers \\
                --project-id my-gcp-project \\
                --dataset my_dataset \\
                --extract-table-prefix my_prefix \\
                --start-bin 0 \\
                --end-bin 100 \\
                --output-bucket-name my_output_bucket \\
                --extract-bucket-path path/to/extract \\
                --obs-columns-to-include "c.cell_type,c.donor_id,c.sex,i.dataset_id"
    """
    bq_ops.extract_bins_in_parallel_workers(
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
