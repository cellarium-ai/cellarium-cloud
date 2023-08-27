"""
Prepare data for extract by randomizing, preprocessing and staging in temporary tables
"""
import argparse
import time
import typing as t

from google.cloud import bigquery

from casp.bq_scripts import constants


def execute_query(client, sql):
    """
    Runs the supplied query
    """
    start = time.time()

    # print(f"RUNNING:\n {sql}")
    query = client.query(sql)
    results = query.result()

    job = client.get_job(query.job_id)
    gb_billed = int(0 if job.total_bytes_billed is None else job.total_bytes_billed) / (1024 * 1024 * 1024)
    print(f"COMPLETED ({time.time() - start} seconds, {gb_billed} GBs scanned)")

    return results


def __get_filter_by_organism_feature_summary_join(
    project: str, dataset: str, filter_by_organism: t.Optional[str]
) -> str:
    if filter_by_organism is None:
        return f"JOIN `{project}.{dataset}.cas_cell_info` c ON (m.cas_cell_index = c.cas_cell_index)"

    return f"""
        JOIN `{project}.{dataset}.cas_cell_info` c ON (m.cas_cell_index = c.cas_cell_index)
            AND c.organism = '{filter_by_organism}'
    """


def __get_filter_by_dataset_feature_summary_join(
    project: str, dataset: str, filter_by_datasets: t.Optional[t.List[str]]
) -> str:
    if filter_by_datasets is None:
        return ""

    filter_by_datasets_str = ", ".join([f"'{x}'" for x in filter_by_datasets])
    return f"""
        JOIN `{project}.{dataset}.cas_ingest_info` i ON (c.cas_ingest_id = i.cas_ingest_id)
             AND i.dataset_filename IN ({filter_by_datasets_str})
    """


def prepare_feature_summary(
    client,
    project,
    dataset,
    extract_table_prefix,
    filter_by_organism: t.Optional[str],
    filter_by_datasets: t.Optional[t.List[str]],
):
    """
    create feature/gene level summary -- scan of full dataset
    """
    filter_by_organism_join = __get_filter_by_organism_feature_summary_join(
        project=project, dataset=dataset, filter_by_organism=filter_by_organism
    )
    filter_by_datasets_join = __get_filter_by_dataset_feature_summary_join(
        project=project, dataset=dataset, filter_by_datasets=filter_by_datasets
    )

    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_feature_summary`
        AS
        SELECT  f.original_feature_id,
                SUM(m.raw_counts) total_raw_counts,
                COUNT(distinct CASE WHEN m.raw_counts > 0 THEN m.cas_cell_index ELSE null END) cells_with_counts
        FROM `{project}.{dataset}.cas_raw_count_matrix` m
        JOIN `{project}.{dataset}.cas_feature_info` f ON (m.cas_feature_index = f.cas_feature_index)
        {filter_by_organism_join}{filter_by_datasets_join}
        GROUP BY f.original_feature_id
    """
    print("Creating Feature Summary...")
    query = execute_query(client, sql)
    return query


def prepare_feature_info(client, project, dataset, extract_table_prefix, fq_allowed_original_feature_ids):
    """
    create subset of features based on
      - a minimum number of observed cells
      - original_feature_id being present in the fully qualified single-column table fq_allowed_original_feature_ids
    and generating with a new feature index value
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_feature_info`
        AS
        SELECT  DENSE_RANK() OVER (ORDER BY fs.feature_name ASC) AS cas_feature_index,
                fs.feature_name as original_feature_id,
                fs.index,
        FROM	`{fq_allowed_original_feature_ids}` fs
        ORDER BY fs.index
    """

    print("Creating Feature Info...")
    query = execute_query(client, sql)
    return query


def __get_prepare_cell_info_join(project: str, dataset: str, filter_by_datasets: t.Optional[t.List[str]]):
    if filter_by_datasets is None:
        return ""

    return f"JOIN `{project}.{dataset}.cas_ingest_info` i ON (i.cas_ingest_id = c.cas_ingest_id)"


def __get_prepare_cell_info_filter(
    datasets: t.Optional[t.List[str]],
    organism: t.Optional[str],
    is_primary_data: t.Optional[bool],
    diseases: t.Optional[t.List[str]],
) -> str:
    if datasets is None and organism is None and is_primary_data is None:
        return ""

    # Parse lists to SQL friendly string
    datasets = ", ".join(f"'{s}'" for s in datasets) if datasets is not None else None
    diseases = ", ".join(f"'{s}'" for s in diseases) if diseases is not None else None
    # Create filters
    filter_by_datasets = f"i.dataset_filename IN ({datasets})" if datasets is not None else None
    filter_by_organism = f"c.organism = '{organism}'" if organism is not None else None
    filter_by_diseases = f"c.disease IN ({diseases})" if diseases is not None else None
    filter_by_is_primary_data = (
        f"c.is_primary_data = {str(is_primary_data).upper()}" if is_primary_data is not None else None
    )
    filters = [filter_by_organism, filter_by_datasets, filter_by_is_primary_data, filter_by_diseases]
    filters = [x for x in filters if x is not None]

    where_body = "\nAND ".join(filters)

    return f"WHERE {where_body}" if where_body != "" else ""


def prepare_cell_info(
    client: "bigquery.Client",
    project: str,
    dataset: str,
    extract_table_prefix: str,
    extract_bin_size: int,
    random_seed_offset: int = 0,
    partition_bin_count: int = 40000,
    partition_size: int = 10,
    filter_by_organism: t.Optional[str] = None,
    filter_by_datasets: t.Optional[t.List[str]] = None,
    filter_by_is_primary_data: t.Optional[bool] = None,
    filter_by_diseases: t.Optional[t.List[str]] = None,
    obs_columns_to_include: t.Optional[t.List[str]] = None,
):
    """
    Randomize cells using farm_fingerprint with an offset so we can have deterministic randomization.  See
    https://towardsdatascience.com/advanced-random-sampling-in-bigquery-sql-7d4483b580bb

    Then allocate cells into bins of `extract_bin_size` cells.  The last bin may be much
    smaller than this requested size, as it is the remainder cells

    :param client: An instance of BigQuery's client object used to interact with the BigQuery service.
    :param project: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param extract_table_prefix: A prefix string for naming tables or for similar purposes.
    :param extract_bin_size: The size for the bins where cells are allocated.
    :param random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization.
        `Defaults:` ``0``
    :param partition_bin_count: The count of bins for partitioning.
        `Default:` ``40000``
    :param partition_size: The size for the partitions.
        `Default:` ``10``
    :param filter_by_organism: Optional filter to specify an organism. If not provided, no filtering occurs.
    :param filter_by_datasets: Optional filter to specify datasets by their names. If not provided, no filtering occurs.
    :param filter_by_is_primary_data: Optional filter to determine if data is primary or not.
    :param filter_by_diseases: Optional filter to specify diseases. If not provided not filtering occurs.
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``. If not
        provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
    """

    join_clause = __get_prepare_cell_info_join(project=project, dataset=dataset, filter_by_datasets=filter_by_datasets)
    where_clause = __get_prepare_cell_info_filter(
        organism=filter_by_organism,
        datasets=filter_by_datasets,
        is_primary_data=filter_by_is_primary_data,
        diseases=filter_by_diseases,
    )
    cas_cell_info_columns = ", ".join([*constants.CAS_CELL_INFO_REQUIRED_COLUMNS, *obs_columns_to_include])
    sql_random_ordering = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_cell_info_randomized`
        AS
        SELECT  {cas_cell_info_columns}
        FROM `{project}.{dataset}.cas_cell_info` c
        {join_clause}
        {where_clause}
        ORDER BY farm_fingerprint(cast(cas_cell_index + {random_seed_offset} as STRING))
    """
    print("Randomizing order of the cells...")
    execute_query(client, sql_random_ordering)

    sql_prepare_cell_info = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_cell_info`
        PARTITION BY RANGE_BUCKET(extract_bin, GENERATE_ARRAY(0,{partition_bin_count},{partition_size}))
        CLUSTER BY extract_bin
        AS
        SELECT  {cas_cell_info_columns},
                CAST(FLOOR((ROW_NUMBER() OVER () - 1) / {extract_bin_size}) as INT) as extract_bin
        FROM `{project}.{dataset}.{extract_table_prefix}__extract_cell_info_randomized` c
    """
    print("Creating Cell Info into extract bins...")
    main_query = execute_query(client, sql_prepare_cell_info)

    sql_remove_random_ordering = f"""
        DROP TABLE `{project}.{dataset}.{extract_table_prefix}__extract_cell_info_randomized`
    """
    print("Removing intermediate table used for random ordering...")
    execute_query(client, sql_remove_random_ordering)
    return main_query


def prepare_extract_matrix(
    client, project, dataset, extract_table_prefix, partition_bin_count=40000, partition_size=10
):
    """
    Create extract table of count data -- remapping feature identifiers, and including batch identifier

    :param partition_bin_count: Number of partition bins that is used by BigQuery to optimize future queries
    :param partition_size: A size of each of the partition bins

    For more info on BigQuery partitioning refer to
    https://cloud.google.com/bigquery/docs/partitioned-tables
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_raw_count_matrix`
        PARTITION BY RANGE_BUCKET(extract_bin, GENERATE_ARRAY(0,{partition_bin_count},{partition_size}))
        CLUSTER BY extract_bin
        AS
        SELECT  b.extract_bin,
                m.cas_cell_index,
                ARRAY_AGG(STRUCT<feature_index int64, raw_counts int64>(ef.cas_feature_index, m.raw_counts)) as feature_data
        FROM `{project}.{dataset}.cas_raw_count_matrix` m
        JOIN `{project}.{dataset}.cas_feature_info` fi ON (m.cas_feature_index = fi.cas_feature_index)
        JOIN `{project}.{dataset}.{extract_table_prefix}__extract_feature_info` ef ON (fi.original_feature_id = ef.original_feature_id)
        JOIN `{project}.{dataset}.{extract_table_prefix}__extract_cell_info` b ON (m.cas_cell_index = b.cas_cell_index)
        GROUP BY 1,2
    """

    print("Creating extract matrix...")
    query = execute_query(client, sql)
    return query


def prepare_extract(
    project: str,
    dataset: str,
    extract_table_prefix: str,
    fq_allowed_original_feature_ids: str,
    extract_bin_size: int,
    ci_random_seed_offset: int = 0,
    ci_partition_bin_count: int = 40000,
    ci_partition_size: int = 10,
    credentials=None,
    filter_by_organism: t.Optional[str] = None,
    filter_by_datasets: t.Optional[t.List[str]] = None,
    filter_by_is_primary_data: t.Optional[bool] = None,
    filter_by_diseases: t.Optional[t.List[str]] = None,
    obs_columns_to_include: t.Optional[t.List[str]] = None,
):
    """
    Prepare CAS BigQuery tables used for data extraction.

    :param project: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param extract_table_prefix: A prefix string for naming tables or for similar purposes.
    :param fq_allowed_original_feature_ids: BigQuery table with feature schema needed for the extract
    :param extract_bin_size: The size for the bins where cells are allocated.
    :param credentials: Google Cloud Service account credentials
    :param ci_random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization
        Used in cas_cell_info table.
        `Defaults:` ``0``
    :param ci_partition_bin_count: The count of bins for partitioning cas_cell_info table.
        `Default:` ``40000``
    :param ci_partition_size: The size for the partitions in cas_cell_info table.
        `Default:` ``10``
    :param filter_by_organism: Optional filter to specify an organism. If not provided, no filtering occurs.
        `Default:` ``None``
    :param filter_by_datasets: Optional filter to specify datasets by their names. If ``None``, no filtering occurs.
        `Default:` ``None``
    :param filter_by_is_primary_data: Optional filter to determine if data is primary or not.
        `Default:` ``None``
    :param filter_by_diseases: Optional filter to specify diseases. If not provided not filtering occurs.
        `Default:` ``None``
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``. If not
        ``None``, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        `Default:` ``None``
    """
    if credentials is None:
        client = bigquery.Client(project=project)
    else:
        client = bigquery.Client(project=project, credentials=credentials)

    prepare_feature_summary(
        client,
        project,
        dataset,
        extract_table_prefix,
        filter_by_organism=filter_by_organism,
        filter_by_datasets=filter_by_datasets,
    )
    prepare_feature_info(client, project, dataset, extract_table_prefix, fq_allowed_original_feature_ids)
    prepare_cell_info(
        client,
        project,
        dataset,
        extract_table_prefix,
        extract_bin_size,
        random_seed_offset=ci_random_seed_offset,
        partition_bin_count=ci_partition_bin_count,
        partition_size=ci_partition_size,
        filter_by_organism=filter_by_organism,
        filter_by_datasets=filter_by_datasets,
        filter_by_is_primary_data=filter_by_is_primary_data,
        filter_by_diseases=filter_by_diseases,
        obs_columns_to_include=obs_columns_to_include,
    )
    prepare_extract_matrix(client, project, dataset, extract_table_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Prepare CASP tables ML Training/Inference Extract"
    )
    parser.add_argument("--project", type=str, help="BigQuery Project", required=True)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix for extract tables", required=True)
    parser.add_argument(
        "--min_observed_cells", type=int, help="minimum observed cells per gene", default=3, required=False
    )
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

    args = parser.parse_args()
    prepare_extract(
        args.project,
        args.dataset,
        args.extract_table_prefix,
        args.min_observed_cells,
        args.fq_allowed_original_feature_ids,
        args.extract_bin_size,
    )
