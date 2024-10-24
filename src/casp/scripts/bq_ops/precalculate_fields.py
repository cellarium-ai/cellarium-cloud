import typing as t

from google.cloud import bigquery

from casp.scripts.bq_ops import constants

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


def precalculate_fields(dataset: str, fields: t.List[str], project: str):
    if not set(fields).issubset(set(SQL_FIELD_MAPPING.keys())):
        raise ValueError("Fields should be one of the `SQL_FIELD_MAPPING` keys")

    client = bigquery.Client(project=project)

    for field in fields:
        print(f"Executing calculation for {field}")
        sql_format = SQL_FIELD_MAPPING[field]
        query_job = client.query(query=sql_format.format(project=project, dataset=dataset))
        _ = query_job.result()  # Wait till command runs
        print("Done.")
