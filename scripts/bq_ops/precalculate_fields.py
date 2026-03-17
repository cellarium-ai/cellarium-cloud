import typing as t

from google.cloud import bigquery

from casp.data_manager import sql
from casp.scripts.bq_ops import constants

PRECALCULATED_FIELD_TEMPLATE_PATH_MAPPING = {
    "total_mrna_umis": constants.TOTAL_MRNA_UMIS_SQL_DIR,
    "bq_row_number": constants.BQ_ROW_NUMBER_SQL_DIR,
}


def precalculate_fields(dataset: str, fields: t.List[str], project: str):
    if not set(fields).issubset(set(PRECALCULATED_FIELD_TEMPLATE_PATH_MAPPING.keys())):
        raise ValueError("Fields should be one of the `SQL_FIELD_MAPPING` keys")

    client = bigquery.Client(project=project)

    for field in fields:
        print(f"Executing calculation for {field}")
        sql_template_path = PRECALCULATED_FIELD_TEMPLATE_PATH_MAPPING[field]

        template_data_rand_ordering = sql.TemplateData(project=project, dataset=dataset)
        sql_query = sql.render(
            template_path=sql_template_path,
            template_data=template_data_rand_ordering,
        )
        query_job = client.query(query=sql_query)
        _ = query_job.result()  # Wait till command runs
        print("Done.")
