from google.cloud import bigquery

MATCH_CELL_RESULTS_SCHEMA = [
    bigquery.SchemaField("query_id", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("match_cas_cell_index", "INT64", mode="REQUIRED"),
    bigquery.SchemaField("match_score", "FLOAT", mode="REQUIRED"),
]
