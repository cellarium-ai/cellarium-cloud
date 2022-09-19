"""
Prepare data for extract by randomizing, preprocessing and staging in temporary tables
"""
import argparse
import time

from google.cloud import bigquery


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


def prepare_feature_summary(client, project, dataset, extract_table_prefix):
    """
    create feature/gene level summary -- scan of full dataset
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_feature_summary`
        AS
        SELECT  f.original_feature_id,
                SUM(m.raw_counts) total_raw_counts,
                COUNT(distinct CASE WHEN m.raw_counts > 0 THEN m.cas_cell_index ELSE null END) cells_with_counts
        FROM `{project}.{dataset}.cas_raw_count_matrix` m
        JOIN `{project}.{dataset}.cas_feature_info` f ON (m.cas_feature_index = f.cas_feature_index)
        GROUP BY f.original_feature_id
    """
    print("Creating Feature Summary...")
    query = execute_query(client, sql)
    return query


def prepare_feature_info(
    client, project, dataset, extract_table_prefix, min_observed_cells, fq_allowed_original_feature_ids
):
    """
    create subset of features based on
      - a minimum number of observed cells
      - original_feature_id being present in the fully qualified single-column table fq_allowed_original_feature_ids
    and generating with a new feature index value
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_feature_info`
        AS
        SELECT  DENSE_RANK() OVER (ORDER BY s.original_feature_id ASC) AS cas_feature_index,
                s.original_feature_id as original_feature_id,
        FROM	`{project}.{dataset}.{extract_table_prefix}__extract_feature_summary` s
        WHERE s.cells_with_counts >= {min_observed_cells}
        AND s.original_feature_id IN (
            SELECT * FROM `{fq_allowed_original_feature_ids}`
        )
        ORDER BY s.original_feature_id
    """

    print("Creating Feature Info...")
    query = execute_query(client, sql)
    return query


def prepare_cell_info(client, project, dataset, extract_table_prefix, extract_bin_size):
    """
    Randomize cells into bins of approximately `extract_bin_size` cells.  The last bin may be much
    smaller than this requested size, as it is the remainder cells
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_cell_info`
        AS
        SELECT  cas_cell_index,
                CAST(FLOOR(rand() * (select count(1) from `{project}.{dataset}.cas_cell_info`) / {extract_bin_size}) as INT) as extract_bin
        FROM	`{project}.{dataset}.cas_cell_info` c
        ORDER BY cas_cell_index
    """

    print("Creating Cell Info and randomizing into extract bins...")
    query = execute_query(client, sql)
    return query


def prepare_extract_matrix(client, project, dataset, extract_table_prefix):
    """
    create extract table of count data -- remapping feature identifiers, and including batch identifier
    """
    sql = f"""
        CREATE OR REPLACE TABLE `{project}.{dataset}.{extract_table_prefix}__extract_raw_count_matrix`
        PARTITION BY RANGE_BUCKET(extract_bin, GENERATE_ARRAY(0,4000,1))
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

    print("Creating Cell Info and randomizing into extract bins...")
    query = execute_query(client, sql)
    return query


def prepare_extract(
    project, dataset, extract_table_prefix, min_observed_cells, fq_allowed_original_feature_ids, extract_bin_size
):
    client = bigquery.Client(project=project)

    prepare_feature_summary(client, project, dataset, extract_table_prefix)
    prepare_feature_info(
        client, project, dataset, extract_table_prefix, min_observed_cells, fq_allowed_original_feature_ids
    )
    prepare_cell_info(client, project, dataset, extract_table_prefix, extract_bin_size)
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