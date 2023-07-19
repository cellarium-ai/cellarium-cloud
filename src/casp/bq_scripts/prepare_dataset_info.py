import numpy as np
import pandas as pd
from google.cloud import bigquery


def prepare_expressed_genes_info(
    project: str, dataset: str, fq_allowed_original_feature_ids: str, credentials=None
) -> pd.DataFrame:
    if credentials is not None:
        client = bigquery.Client(project=project, credentials=credentials)
    else:
        client = bigquery.Client(project=project)
    sql_expressed_genes = f"""
        SELECT  cas_ingest_id,
                ARRAY_AGG(original_feature_id) as expressed_genes
        FROM `{project}.{dataset}.cas_feature_info`
        GROUP BY 1
    """
    sql_get_all_features = f"""
        SELECT feature_name FROM `{fq_allowed_original_feature_ids}`
    """
    sql_get_number_of_ingests = f"""
        SELECT COUNT(*) AS number_of_ingests FROM `{project}.{dataset}.cas_ingest_info`
    """
    expressed_genes_result = client.query(sql_expressed_genes)
    all_features_result = client.query(sql_get_all_features)
    number_of_ingests_result = client.query(sql_get_number_of_ingests)

    all_features = np.array([x["feature_name"] for x in all_features_result])
    number_of_ingests = next(iter(number_of_ingests_result))["number_of_ingests"]

    gene_expression_existance_mask = np.zeros((number_of_ingests, len(all_features)))
    cas_ingest_ids = []

    for i, row in enumerate(expressed_genes_result):
        genes_expressed_in_ingest = row["expressed_genes"]
        cas_ingest_ids.append(row["cas_ingest_id"])
        gene_expression_existance_mask[i, :] = np.isin(all_features, genes_expressed_in_ingest)

    return pd.DataFrame(
        data=gene_expression_existance_mask, index=pd.Index(cas_ingest_ids, name="cas_ingest_id"), columns=all_features
    )


def prepare_all_cell_types(project: str, dataset: str, credentials) -> pd.DataFrame:
    if credentials is not None:
        client = bigquery.Client(project=project, credentials=credentials)
    else:
        client = bigquery.Client(project=project)

    sql_get_all_cel_types = f"""
        SELECT DISTINCT(cell_type)
        FROM `{project}.{dataset}.cas_cell_info`
    """
    return client.query(sql_get_all_cel_types).to_dataframe()
