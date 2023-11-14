import typing as t

from google.cloud import bigquery

from casp.datastore_manager import sql
from casp.services import settings, utils
from casp.services.api import constants, schemas
from casp.services.db import models

SQL_GET_ALL_GENE_SCHEMAS = """
    SELECT table_name FROM `{project}.cas_reference_data.INFORMATION_SCHEMA.TABLES`;
"""

SQL_GET_SCHEMA_BY_NAME = """
    SELECT feature_name
    FROM `{project}.cas_reference_data.{schema_name}`
    ORDER BY index;
"""


def get_application_info() -> schemas.ApplicationInfo:
    """
    :return: Object with CAS application information
    """
    return schemas.ApplicationInfo(
        default_feature_schema=settings.DEFAULT_FEATURE_SCHEMA, application_version=settings.APP_VERSION
    )


def get_feature_schemas() -> t.List[schemas.FeatureSchemaInfo]:
    """
    :return: List of gene schema objects
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    sql_query = SQL_GET_ALL_GENE_SCHEMAS.format(project=project)
    query_result = client.query(sql_query).result()
    return [schemas.FeatureSchemaInfo(schema_name=row["table_name"]) for row in query_result]


def get_feature_schema_by(schema_name: str) -> t.List[str]:
    """
    Get ensebl gene ids that are in schema
    :param schema_name: Unique schema name
    :return: List of Ensebl gene ids
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    sql_query = SQL_GET_SCHEMA_BY_NAME.format(project=project, schema_name=schema_name)
    query_result = client.query(sql_query).result()
    return [row["feature_name"] for row in query_result]


def _parse_match_object(row) -> t.Dict[str, t.Any]:
    """
    Parse a single row of a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param row: Row of the BigQuery job object representing the query execution.
    :return: Dictionary representing the query results.
    """
    return {
        "cell_type": row["cell_type"],
        "cell_count": row["cell_count"],
        "min_distance": row["min_distance"],
        "p25_distance": row["p25_distance"],
        "median_distance": row["median_distance"],
        "p75_distance": row["p75_distance"],
        "max_distance": row["max_distance"],
    }


def _parse_match_object_dev_details(row) -> t.Dict[str, t.Any]:
    """
    Parse a single row of a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param row: Row of the BigQuery job object representing the query execution.
    :return: Dictionary representing the query results.
    """
    dataset_ids_with_counts = []
    for dataset_ids_with_counts_struct in row["dataset_ids_with_counts"]:
        dataset_ids_with_counts.append(
            {
                "dataset_id": dataset_ids_with_counts_struct["dataset_id"],
                "count_per_dataset": dataset_ids_with_counts_struct["count_per_dataset"],
                "min_distance": dataset_ids_with_counts_struct["min_distance"],
                "max_distance": dataset_ids_with_counts_struct["max_distance"],
                "median_distance": dataset_ids_with_counts_struct["median_distance"],
                "mean_distance": dataset_ids_with_counts_struct["mean_distance"],
            }
        )
    return {
        "cell_type": row["cell_type"],
        "cell_count": row["cell_count"],
        "min_distance": row["min_distance"],
        "p25_distance": row["p25_distance"],
        "median_distance": row["median_distance"],
        "p75_distance": row["p75_distance"],
        "max_distance": row["max_distance"],
        "dataset_ids_with_counts": dataset_ids_with_counts,
    }


def _parse_match_query_job(query_job, include_dev_details: bool = False) -> t.List[t.Dict[str, t.Any]]:
    """
    Parses a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param query_job: BigQuery job object representing the query execution.
    :param include_dev_details: Boolean indicating whether to include a breakdown of the number of cells by dataset
    :return: List of dictionaries representing the query results.
    """
    results = []

    last_query_id = None
    data = {}

    for row in query_job:
        if last_query_id is None or last_query_id != row["query_id"]:
            # emit data and reset state if this isn't the first time through
            if last_query_id is not None:
                results.append(data)

            data = {"query_cell_id": row["query_id"], "matches": []}
            last_query_id = data["query_cell_id"]

        if not include_dev_details:
            x = _parse_match_object(row=row)
        else:
            x = _parse_match_object_dev_details(row=row)

        data["matches"].append(x)

    results.append(data)
    return results


def get_match_query_metadata(cas_model: models.CASModel, match_temp_table_fqn: str) -> t.List[t.Dict[str, t.Any]]:
    """
    Executes a BigQuery query to retrieve metadata for a matching query.

    :param cas_model: The CASModel containing dataset information
    :param match_temp_table_fqn: The fully-qualified name of the temporary table.
    :return: The BigQuery job object representing the query execution.
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    sql_template_data = sql.TemplateData(
        project=project, dataset=cas_model.bq_dataset_name, temp_table_fqn=match_temp_table_fqn
    )
    sql_query = sql.render(constants.MATCH_METADATA_SQL_TMPL_PATH, template_data=sql_template_data)
    query_job = client.query(sql_query)
    return _parse_match_query_job(query_job=query_job)


def get_match_query_metadata_dev_details(
    cas_model: models.CASModel, match_temp_table_fqn: str
) -> t.List[t.Dict[str, t.Any]]:
    """
    Executes a BigQuery query to retrieve metadata for a matching query. The returned query, similar to
    :func:`get_match_query_metadata`, includes a breakdown of the number of cells that matched each cell type by dataset.

    :param cas_model: The CASModel containing dataset information
    :param match_temp_table_fqn: The fully-qualified name of the temporary table.
    :return: The BigQuery job object representing the query execution.
    """
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    sql_template_data = sql.TemplateData(
        project=project, dataset=cas_model.bq_dataset_name, temp_table_fqn=match_temp_table_fqn
    )
    sql_query = sql.render(constants.MATCH_METADATA_DEV_DETAILS_SQL_TMPL_PATH, template_data=sql_template_data)
    query_job = client.query(sql_query)
    return _parse_match_query_job(query_job=query_job, include_dev_details=True)
