import typing as t

from google.cloud import bigquery

from casp.services import utils

SQL_GET_ORIGINAL_CELL_IDS = """
    SELECT cas_cell_index, original_cell_id from `{project}.{dataset}.cas_cell_info`
    WHERE cas_cell_index IN ({cas_cell_ids_clause});
"""


def get_original_cell_ids_by(dataset_name: str, cas_cell_ids: t.List[str]) -> t.List[t.Tuple[str, str]]:
    """
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
