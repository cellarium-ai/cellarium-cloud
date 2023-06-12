import typing as t

from google.cloud import bigquery

from casp.services import settings, utils
from casp.services.api import schemas

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
