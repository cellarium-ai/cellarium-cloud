import typing as t

import pandas as pd
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

SQL_GET_CAS_CELL_INFO_COLUMNS = """
    SELECT {columns_to_select} from `{project}.{dataset}.cas_cell_info` c
    JOIN `{project}.{dataset}.cas_ingest_info` i ON (c.cas_ingest_id = i.cas_ingest_id)
    WHERE cas_cell_index IN ({cas_cell_ids_clause});
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


def get_cas_cell_info_columns_df(
    dataset_name: str, cas_cell_ids: t.List[str], columns_to_select: t.List[str]
) -> "pd.DataFrame":
    """
    Retrieve selected columns for specific cell IDs from the `cas_cell_info` table in BigQuery.

    Given a dataset name, a list of cell IDs (`cas_cell_ids`), and a list of columns to select, this function
    formulates an SQL query to fetch the relevant data from the `cas_cell_info` table in BigQuery. The resulting DataFrame
    has its index set to the `cas_cell_index`, ensuring the order matches the input `cas_cell_ids` list.

    :param dataset_name: Name of the dataset in BigQuery to fetch data from.
    :param cas_cell_ids: List of cell IDs for which to fetch data.
    :param columns_to_select: List of column names to be selected from the `cas_cell_info` table.
    :return: :class:`pd.DataFrame` instance containing the selected columns for the specified cell IDs, indexed by
        `cas_cell_index`.
    """
    credentials, project = utils.get_google_service_credentials()
    cas_cell_ids_str_clause = ",".join(cas_cell_ids)
    columns_to_select_str_clause = ", ".join(["cas_cell_index", *columns_to_select])
    sql_query = SQL_GET_CAS_CELL_INFO_COLUMNS.format(
        project=project,
        dataset=dataset_name,
        cas_cell_ids_clause=cas_cell_ids_str_clause,
        columns_to_select=columns_to_select_str_clause,
    )
    df = pd.read_gbq(query=sql_query, credentials=credentials, project_id=project)
    df = df.set_index("cas_cell_index")
    # Casting to string, because extract chunks indexes are string, while BQ client
    # returns integer as it is an integer in the database schema definition
    df.index = df.index.astype(str)

    df.index.names = [""]
    # Order cas_cell_info according to input cas_cell_ids order, because BQ doesn't keep the order
    df = df.loc[cas_cell_ids, :]
    return df
