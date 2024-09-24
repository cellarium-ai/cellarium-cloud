import typing as t

import numpy as np
import pandas as pd
from google.cloud import bigquery
from google.oauth2.service_account import Credentials
from casp.data_manager import sql
from casp.scripts.bq_ops import constants


def prepare_measured_genes_info(
    project: str, dataset: str, fq_allowed_original_feature_ids: str, credentials: t.Optional[Credentials] = None
) -> pd.DataFrame:
    if credentials is not None:
        client = bigquery.Client(project=project, credentials=credentials)
    else:
        client = bigquery.Client(project=project)
    sql_measured_genes = f"""
        SELECT  cas_ingest_id,
                ARRAY_AGG(original_feature_id) as measured_genes
        FROM `{project}.{dataset}.cas_feature_info`
        GROUP BY 1
    """
    sql_get_all_features = f"""
        SELECT feature_name FROM `{fq_allowed_original_feature_ids}`
    """
    sql_get_number_of_ingests = f"""
        SELECT COUNT(*) AS number_of_ingests FROM `{project}.{dataset}.cas_ingest_info`
    """
    measured_genes_result = client.query(sql_measured_genes)
    all_features_result = client.query(sql_get_all_features)
    number_of_ingests_result = client.query(sql_get_number_of_ingests)

    all_features = np.array([x["feature_name"] for x in all_features_result])
    number_of_ingests = next(iter(number_of_ingests_result))["number_of_ingests"]

    gene_expression_existance_mask = np.zeros((number_of_ingests, len(all_features)))
    cas_ingest_ids = []

    for i, row in enumerate(measured_genes_result):
        genes_measured_in_ingest = row["measured_genes"]
        cas_ingest_ids.append(row["cas_ingest_id"])
        gene_expression_existance_mask[i, :] = np.isin(all_features, genes_measured_in_ingest)

    return pd.DataFrame(
        data=gene_expression_existance_mask, index=pd.Index(cas_ingest_ids, name="cas_ingest_id"), columns=all_features
    )


def prepare_categorical_variables(project: str, dataset: str, extract_table_prefix: str) -> t.Dict[str, t.List[str]]:
    client = bigquery.Client(project=project)

    table_ref = bigquery.TableReference.from_string(f"{project}.{dataset}.cas_cell_info")
    table = client.get_table(table_ref)

    columns_unique_values = {}

    for schema_field in table.schema:
        if schema_field.field_type == "STRING":
            template_data = sql.TemplateData(
                project=project,
                dataset=dataset,
                extract_table_prefix=extract_table_prefix,
                column_name=schema_field.name,
            )

            sql_query = sql.render(
                template_path=constants.PREPARE_CATEGORICAL_VARIABLE_SQL_DIR, template_data=template_data
            )

            columns_unique_values[schema_field.name] = [x[schema_field.name] for x in client.query(sql_query).result()]

    return columns_unique_values


# def prepare_all_cell_types(project: str, dataset: str, credentials: t.Optional[Credentials] = None) -> pd.DataFrame:
#     if credentials is not None:
#         client = bigquery.Client(project=project, credentials=credentials)
#     else:
#         client = bigquery.Client(project=project)
#
#     sql_get_all_cel_types = f"""
#         SELECT DISTINCT(cell_type)
#         FROM `{project}.{dataset}.cas_cell_info`
#     """
#     return client.query(sql_get_all_cel_types).to_dataframe()
