import typing as t

from google.cloud import bigquery
from casp.scripts.bq_ops import constants


def get_measured_genes_info_filepath(bucket_name: str, extract_bucket_path: str) -> str:
    return (
        f"gs://{bucket_name}/{extract_bucket_path}"
        f"/{constants.SHARED_META_DIR_NAME}"
        f"/{constants.MEASURED_GENES_INFO_FILE_NAME}"
    )


def get_categorical_columns_metadata_filepath(bucket_name: str, extract_bucket_path: str) -> str:
    return (
        f"gs://{bucket_name}/{extract_bucket_path}"
        f"/{constants.SHARED_META_DIR_NAME}"
        f"/{constants.CATEGORICAL_COLUMNS_META_FILE_NAME}"
    )


def get_categorical_column_names_from_cell_info_table(
    project: str, dataset: str, extract_table_prefix: str
) -> t.List[str]:
    client = bigquery.Client(project=project)

    table_ref = bigquery.TableReference.from_string(f"{project}.{dataset}.{extract_table_prefix}__extract_cell_info")
    table = client.get_table(table_ref)

    return [schema_field.name for schema_field in table.schema if schema_field.field_type == "STRING"]
