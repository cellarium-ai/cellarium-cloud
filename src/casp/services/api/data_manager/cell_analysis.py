import typing as t
import uuid
from datetime import datetime, timedelta

from google.cloud import bigquery
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor

from casp.data_manager import BaseDataManager, sql
from casp.services import settings
from casp.services.api.data_manager import bigquery_response_parsers, bigquery_schemas
from casp.services.db import models


class CellAnalysisDataManager(BaseDataManager):
    """
    Data Manager for accessing data withing the framework of Cell Analysis.
    """

    # Directories for SQL templates
    KNN_QUERY_TEMPLATE_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_templates/knn_query_match"

    # SQL template file paths
    SQL_MATCH_METADATA = f"{KNN_QUERY_TEMPLATE_DIR}/match_metadata.sql.mako"
    SQL_MATCH_METADATA_DEV_DETAILS = f"{KNN_QUERY_TEMPLATE_DIR}/match_metadata_dev_details.sql.mako"

    def insert_matches_to_temp_table(self, query_ids: t.List[str], knn_response: t.List[t.List[MatchNeighbor]]) -> str:
        """
        Insert matches to temporary table in async manner in a separate thread.

        Note: This function executes I/O-bound operation (BigQuery insert) in a separate thread to avoid blocking the
        main thread. Parsing of the response is done in the main thread in sync manner.

        :param query_ids: List of query ids (original cell ids from the input file).
        :param knn_response: List of lists of MatchNeighbor objects returned by space vector search service.

        :return: The fully-qualified name of the temporary table.
        """
        my_uuid = str(uuid.uuid4())[:8]
        temp_table_fqn = f"{settings.API_REQUEST_TEMP_TABLE_DATASET}.api_request_matches_{my_uuid}"
        table = bigquery.Table(temp_table_fqn, schema=bigquery_schemas.MATCH_CELL_RESULTS_SCHEMA)
        table.expires = datetime.now() + timedelta(minutes=settings.API_REQUEST_TEMP_TABLE_DATASET_EXPIRATION)

        self.bigquery_client.create_table(table)

        rows_to_insert = []
        for i in range(0, len(knn_response)):
            query_id = query_ids[i]
            for match in knn_response[i]:
                rows_to_insert.append(
                    {
                        "query_id": str(query_id),
                        "match_cas_cell_index": int(match.id),
                        "match_score": float(match.distance),
                    }
                )

        job_config = bigquery.LoadJobConfig(schema=bigquery_schemas.MATCH_CELL_RESULTS_SCHEMA)

        job = self.bigquery_client.load_table_from_json(
            json_rows=rows_to_insert, destination=temp_table_fqn, job_config=job_config
        )

        job.result()  # Wait for the job to complete

        return temp_table_fqn

    def get_match_query_metadata(
        self, cas_model: models.CASModel, match_temp_table_fqn: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Execute a BigQuery query to retrieve metadata for a matching query.

        :param cas_model: The CASModel containing dataset information
        :param match_temp_table_fqn: The fully-qualified name of the temporary table.

        :return: The BigQuery job object representing the query execution.
        """
        sql_template_data = sql.TemplateData(
            project=self.project, dataset=cas_model.bq_dataset_name, temp_table_fqn=match_temp_table_fqn
        )
        sql_query = sql.render(self.SQL_MATCH_METADATA, template_data=sql_template_data)

        query_job = self.bigquery_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job)

    def get_match_query_metadata_dev_details(
        self, cas_model: models.CASModel, match_temp_table_fqn: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Execute a BigQuery query to retrieve metadata for a matching query. The returned query, similar to
        :meth:`get_match_query_metadata`, includes a breakdown of the number of cells that matched each cell type
        by dataset.

        :param cas_model: The CASModel containing dataset information
        :param match_temp_table_fqn: The fully-qualified name of the temporary table.

        :return: The BigQuery job object representing the query execution.
        """
        sql_template_data = sql.TemplateData(
            project=self.project, dataset=cas_model.bq_dataset_name, temp_table_fqn=match_temp_table_fqn
        )
        sql_query = sql.render(self.SQL_MATCH_METADATA_DEV_DETAILS, template_data=sql_template_data)

        query_job = self.bigquery_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job, include_dev_details=True)
