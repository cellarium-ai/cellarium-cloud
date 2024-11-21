"""
Prepare data for extract by randomizing, preprocessing and staging in temporary tables
"""

import json
import pickle
import time
import typing as t

from google.cloud import bigquery
from google.oauth2.service_account import Credentials
from smart_open import open

from casp.data_manager import sql
from casp.scripts.bq_ops import constants, extract_metadata_utils
from casp.scripts.bq_ops.prepare_dataset_info import prepare_categorical_variables, prepare_measured_genes_info


def execute_query(client: bigquery.Client, sql: str):
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


def prepare_feature_summary(
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    filters: t.Optional[t.Dict[str, t.Any]],
):
    """
    create feature/gene level summary -- scan of full dataset
    """
    template_data = sql.TemplateData(
        project=project, dataset=dataset, filters=filters, extract_table_prefix=extract_table_prefix
    )
    sql_query = sql.render(template_path=constants.PREPARE_FEATURE_SUMMARY_TEMPLATE_DIR, template_data=template_data)
    print("Creating Feature Summary...")
    query = execute_query(client, sql_query)
    return query


def prepare_feature_info(
    client: bigquery.Client, project: str, dataset: str, extract_table_prefix: str, fq_allowed_original_feature_ids: str
):
    """
    create subset of features based on
      - a minimum number of observed cells
      - original_feature_id being present in the fully qualified single-column table fq_allowed_original_feature_ids
    and generating with a new feature index value
    """
    sql_query = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_feature_info`
        AS
        SELECT  DENSE_RANK() OVER (ORDER BY fs.feature_name ASC) AS cas_feature_index,
                fs.feature_name as original_feature_id,
                fs.index,
        FROM	`{fq_allowed_original_feature_ids}` fs
        ORDER BY fs.index
    """

    print("Creating Feature Info...")
    query = execute_query(client, sql_query)
    return query


def prepare_cell_info(
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    extract_bin_size: t.Optional[int] = None,
    assign_bin_by_category: bool = False,
    extract_bin_category_column_name: t.Optional[str] = None,
    random_seed_offset: int = 0,
    partition_bin_count: int = 40000,
    partition_size: int = 10,
    filters: t.Optional[t.Dict[str, str]] = None,
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
    :param extract_bin_size: The size for the bins where cells are allocated. Used only if `assign_bin_by_category`
        is False |br|
        `Default:` ``None``
    :param assign_bin_by_category: Whether ``extract_bin`` has to be assigned as a label to a categorical column |br|
        `Default:` ``False``
    :param extract_bin_category_column_name: Which categorical column to use for labeling the bins. Used only if
        `assign_bin_by_category` is True. |br|
        `Default:` ``None``
    :param random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization. |br|
        `Default:` ``0``
    :param partition_bin_count: The count of bins for partitioning. |br|
        `Default:` ``40000``
    :param partition_size: The size for the partitions. |br|
        `Default:` ``10``
    :param filters: Filters that have to be included in a SQL query. Filter format should follow convention described in
        :func:`casp.bq_manager.sql.mako_helpers.parse_where_body`
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    if not assign_bin_by_category and extract_bin_size is None:
        raise ValueError("If `assign_bin_by_category` is False, `extract_bin_size` must be provided.")
    if assign_bin_by_category and extract_bin_category_column_name is None:
        raise ValueError("If `assign_bin_by_category` is True, `extract_bin_category_column_name` must be provided.")

    template_data_rand_ordering = sql.TemplateData(
        project=project,
        dataset=dataset,
        select=obs_columns_to_include,
        filters=filters,
        random_seed_offset=random_seed_offset,
        extract_table_prefix=extract_table_prefix,
    )
    sql_random_ordering = sql.render(
        template_path=constants.PREPARE_CELL_INFO_RAND_TEMPLATE_DIR,
        template_data=template_data_rand_ordering,
    )
    print("Randomizing order of the cells...")
    execute_query(client, sql_random_ordering)

    template_data_prepare_cell_info = sql.TemplateData(
        project=project,
        dataset=dataset,
        select=obs_columns_to_include,
        filters=filters,
        partition_bin_count=partition_bin_count,
        partition_size=partition_size,
        extract_bin_size=extract_bin_size,
        assign_bin_by_category=assign_bin_by_category,
        extract_bin_category_column_name=extract_bin_category_column_name,
        extract_table_prefix=extract_table_prefix,
    )
    sql_prepare_cell_info = sql.render(
        template_path=constants.PREPARE_CELL_INFO_TEMPLATE_DIR,
        template_data=template_data_prepare_cell_info,
    )
    print("Creating Cell Info into extract bins...")
    main_query = execute_query(client, sql_prepare_cell_info)
    template_data_drop = sql.TemplateData(project=project, dataset=dataset, extract_table_prefix=extract_table_prefix)
    sql_remove_random_ordering = sql.render(
        template_path=constants.DROP_PREPARE_CI_RAND_TEMPLATE_DIR, template_data=template_data_drop
    )
    print("Removing intermediate table used for random ordering...")
    execute_query(client, sql_remove_random_ordering)
    return main_query


def prepare_extract_matrix(
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    partition_bin_count: int = 40000,
    partition_size: int = 10,
):
    """
    Create extract table of count data -- remapping feature identifiers, and including batch identifier

    :param client: BigQuery client.
    :param project: Google Cloud project.
    :param dataset: BigQuery dataset.
    :param extract_table_prefix: A prefix string for naming tables for the extract.
    :param partition_bin_count: Number of partition bins that is used by BigQuery to optimize future queries.
    :param partition_size: A size of each of the partition bins.

    For more info on BigQuery partitioning refer to
    https://cloud.google.com/bigquery/docs/partitioned-tables
    """
    sql_query = f"""
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
    query = execute_query(client, sql_query)
    return query


def prepare_extract_tables(
    project: str,
    dataset: str,
    extract_table_prefix: str,
    fq_allowed_original_feature_ids: str,
    extract_bin_size: t.Optional[int] = None,
    assign_bin_by_category: bool = False,
    extract_bin_category_column_name: t.Optional[str] = None,
    ci_random_seed_offset: int = 0,
    ci_partition_bin_count: int = 40000,
    ci_partition_size: int = 10,
    credentials: t.Optional[Credentials] = None,
    filters: t.Optional[t.Dict[str, t.Any]] = None,
    obs_columns_to_include: t.Optional[t.List[str]] = None,
):
    """
    Prepare CAS BigQuery tables used for data extraction.

    :param project: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param extract_table_prefix: A prefix string for naming tables or for similar purposes.
    :param fq_allowed_original_feature_ids: BigQuery table with feature schema needed for the extract
    :param extract_bin_size: The size for the bins where cells are allocated. Used only if `assign_bin_by_category`
        is False |br|
        `Default:` ``None``
    :param assign_bin_by_category: Whether ``extract_bin`` has to be assigned as a label to a categorical column |br|
        `Default:` ``False``
    :param extract_bin_category_column_name: Which categorical column to use for labeling the bins. Used only if
        `assign_bin_by_category` is True. |br|
        `Default:` ``None``
    :param credentials: Google Cloud Service account credentials
    :param ci_random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization
        Used in cas_cell_info table.
        `Defaults:` ``0``
    :param ci_partition_bin_count: The count of bins for partitioning cas_cell_info table.
        `Default:` ``40000``
    :param ci_partition_size: The size for the partitions in cas_cell_info table.
        `Default:` ``10``
    :param filters: Filters that have to be included in a SQL query. Filter format should follow convention described in
        :func:`casp.bq_manager.sql.mako_helpers.parse_where_body`
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    if credentials is None:
        client = bigquery.Client(project=project)
    else:
        client = bigquery.Client(project=project, credentials=credentials)

    prepare_feature_summary(client, project, dataset, extract_table_prefix, filters=filters)
    prepare_feature_info(client, project, dataset, extract_table_prefix, fq_allowed_original_feature_ids)
    prepare_cell_info(
        client=client,
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        extract_bin_size=extract_bin_size,
        assign_bin_by_category=assign_bin_by_category,
        extract_bin_category_column_name=extract_bin_category_column_name,
        random_seed_offset=ci_random_seed_offset,
        partition_bin_count=ci_partition_bin_count,
        partition_size=ci_partition_size,
        filters=filters,
        obs_columns_to_include=obs_columns_to_include,
    )
    prepare_extract_matrix(client, project, dataset, extract_table_prefix)


def prepare_extract(
    project_id: str,
    dataset: str,
    extract_table_prefix: str,
    fq_allowed_original_feature_ids: str,
    filters_json_path: str,
    obs_columns_to_include: str,
    bucket_name: str,
    extract_bucket_path: str,
    extract_bin_size: t.Optional[int] = None,
    assign_bin_by_category: bool = False,
    extract_bin_category_column_name: t.Optional[str] = None,
    ci_random_seed_offset: int = 0,
    ci_partition_bin_count: int = 40000,
    ci_partition_size: int = 10,
):
    """
    Prepare CAS BigQuery tables used for data extraction.

    :param project_id: The ID of the Google Cloud project where the BigQuery dataset is hosted.
    :param dataset: The ID of the dataset in BigQuery where the data is hosted.
    :param extract_table_prefix: A prefix string for naming tables or for similar purposes.
    :param fq_allowed_original_feature_ids: BigQuery table with feature schema needed for the extract
    :param bucket_name: GCS Bucket name where to store the metadata files.
    :param extract_bucket_path: GCS Bucket path where the extract will be executed. Used to save metadata files there.
        It is required to use the same bucket path during extract
    :param extract_bin_size: The size for the bins where cells are allocated. Used only if `assign_bin_by_category`
        is False |br|
        `Default:` ``None``
    :param assign_bin_by_category: Whether ``extract_bin`` has to be assigned as a label to a categorical column |br|
        `Default:` ``False``
    :param extract_bin_category_column_name: Which categorical column to use for labeling the bins. Used only if
        `assign_bin_by_category` is True. |br|
        `Default:` ``None``
    :param ci_random_seed_offset: Optional offset for the farm_fingerprint for deterministic randomization
        Used in cas_cell_info table.
        `Defaults:` ``0``
    :param ci_partition_bin_count: The count of bins for partitioning cas_cell_info table.
        `Default:` ``40000``
    :param ci_partition_size: The size for the partitions in cas_cell_info table.
        `Default:` ``10``
    :param filters_json_path: A path to json with filters. Filters that have to be included in a SQL query.
    :param obs_columns_to_include: Comma separated string of columns from `cas_cell_info` table to include in
        ``adata.obs``. If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    with open(filters_json_path) as f:
        filters = json.loads(f.read())

    obs_columns_to_include = obs_columns_to_include.split(",")

    prepare_extract_tables(
        project=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        assign_bin_by_category=assign_bin_by_category,
        extract_bin_category_column_name=extract_bin_category_column_name,
        ci_random_seed_offset=ci_random_seed_offset,
        ci_partition_bin_count=ci_partition_bin_count,
        ci_partition_size=ci_partition_size,
        filters=filters,
        obs_columns_to_include=obs_columns_to_include,
    )
    print("Preparing measured genes info...")
    measured_genes_info_df = prepare_measured_genes_info(
        project=project_id,
        dataset=dataset,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
    )
    print("Preparing categorical columns info...")
    categorical_columns_metadata = prepare_categorical_variables(
        project=project_id, dataset=dataset, extract_table_prefix=extract_table_prefix
    )
    measured_genes_info_filepath = extract_metadata_utils.get_measured_genes_info_filepath(
        bucket_name=bucket_name, extract_bucket_path=extract_bucket_path
    )
    categorical_columns_meta_filepath = extract_metadata_utils.get_categorical_columns_metadata_filepath(
        bucket_name=bucket_name, extract_bucket_path=extract_bucket_path
    )

    print("Uploading measured genes info...")
    measured_genes_info_df.to_csv(path_or_buf=measured_genes_info_filepath)

    print("Uploading categorical columns info...")
    with open(categorical_columns_meta_filepath, "wb") as f:
        pickle.dump(obj=categorical_columns_metadata, file=f)
