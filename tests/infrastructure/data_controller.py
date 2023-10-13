import typing as t

from google.cloud import bigquery

from casp.services import utils

SQL_GET_ORIGINAL_CELL_IDS = """
    SELECT cas_cell_index, original_cell_id from `{project}.{dataset}.cas_cell_info`
    WHERE cas_cell_index IN ({cas_cell_ids_clause});
"""

SQL_GET_RAW_MRNA_COUNTS = """
    SELECT m.cas_cell_index, SUM(raw_counts) AS raw_mrna_counts
    FROM `{project}.{dataset}.cas_raw_count_matrix` m
    JOIN `{project}.{dataset}.cas_feature_info` fi ON (fi.cas_feature_index = m.cas_feature_index)
    JOIN `{project}.{dataset}.cas_cell_info` ci ON (ci.cas_cell_index = m.cas_cell_index)
        WHERE fi.feature_biotype = 'gene'
        AND ci.cas_cell_index IN ({cas_cell_ids_clause})
    GROUP BY cas_cell_index;
"""


def get_original_cell_ids_by(dataset_name: str, cas_cell_ids: t.List[str]) -> t.List[t.Tuple[str, str]]:
    """
    Get original cell ids by CAS cell indexes

    :param dataset_name: Name of the BigQuery Dataset
    :param cas_cell_ids: List of cas cell ids to filter by
    :return: List of original cell ids form the initial data source mapped from the BQ
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    cas_cell_ids_str_clause = ",".join(cas_cell_ids)
    sql_query = SQL_GET_ORIGINAL_CELL_IDS.format(
        project=project, dataset=dataset_name, cas_cell_ids_clause=cas_cell_ids_str_clause
    )
    query_result = client.query(sql_query).result()

    return [(str(row["cas_cell_index"]), row["original_cell_id"]) for row in query_result]


def get_raw_mrna_counts_by(dataset_name: str, cas_cell_ids: t.List[str]) -> t.List[t.Tuple[str, int]]:
    """
    Get mRNA counts from BigQuery by CAS cell indexes

    :param dataset_name: Name of the BigQuery Dataset
    :param cas_cell_ids: List of cas cell ids to filter by
    :return: List of tuples, each containing (cas_cell_index, calculated_mrna_counts).
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    cas_cell_ids_str_clause = ",".join(cas_cell_ids)

    sql_query = SQL_GET_RAW_MRNA_COUNTS.format(
        project=project, dataset=dataset_name, cas_cell_ids_clause=cas_cell_ids_str_clause
    )
    query_result = client.query(sql_query).result()

    return [(str(row["cas_cell_index"]), row["raw_mrna_counts"]) for row in query_result]
