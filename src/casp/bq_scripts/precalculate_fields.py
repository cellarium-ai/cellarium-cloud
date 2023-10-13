import typing as t

from google.cloud import bigquery

from casp.bq_scripts import constants
from casp.services.utils import get_google_service_credentials

SQL_PRECALCULATE_TOTAL_MRNA_UMIS_FORMAT = f"""
    UPDATE `{{project}}.{{dataset}}.cas_cell_info` ci
    SET ci.total_mrna_umis = m.raw_counts_total
    FROM (
        SELECT cas_cell_index, SUM(raw_counts) AS raw_counts_total
        FROM `{{project}}.{{dataset}}.cas_raw_count_matrix` m
        JOIN `{{project}}.{{dataset}}.cas_feature_info` fi ON (fi.cas_feature_index = m.cas_feature_index)
        WHERE fi.feature_biotype = '{constants.MRNA_FEATURE_BIOTYPE_NAME}'
        GROUP BY cas_cell_index
    ) AS m
    WHERE ci.cas_cell_index = m.cas_cell_index;
"""

SQL_FIELD_MAPPING = {"total_mrna_umis": SQL_PRECALCULATE_TOTAL_MRNA_UMIS_FORMAT}


def precalculate_fields(dataset: str, fields: t.List[str]):
    assert set(fields).issubset(set(SQL_FIELD_MAPPING.keys())), "Fields should be one of `SQL_FIELD_MAPPING` key"
    credentials, project = get_google_service_credentials()
    client = bigquery.Client(credentials=credentials)

    for field in fields:
        sql_format = SQL_FIELD_MAPPING[field]
        client.query(query=sql_format.format(project=project, dataset=dataset))
