import typing as t
import uuid
from datetime import datetime, timedelta

import sqlalchemy as sa
from google.cloud import bigquery

from casp.data_manager import BaseDataManager, sql
from casp.services import settings
from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.data_manager import bigquery_response_parsers, bigquery_schemas, exceptions
from casp.services.db import models


class CellOperationsDataManager(BaseDataManager):
    """
    Data Manager for making data operations in Cellarium Cloud storage.
    """

    # Directories for SQL templates
    CELL_OPERATIONS_TEMPLATE_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_templates/cell_operations"

    # SQL template file paths
    SQL_MATCH_METADATA = f"{CELL_OPERATIONS_TEMPLATE_DIR}/get_neighborhood_distance_summary.sql.mako"
    SQL_MATCH_METADATA_DEV_DETAILS = (
        f"{CELL_OPERATIONS_TEMPLATE_DIR}/get_neighborhood_distance_summary_dev_details.sql.mako"
    )

    def insert_matches_to_temp_table(self, query_ids: t.List[str], knn_response: MatchResult) -> str:
        """
        Insert matches to temporary table in async manner in a separate thread.

        Note: This function executes I/O-bound operation (BigQuery insert) in a separate thread to avoid blocking the
        main thread. Parsing of the response is done in the main thread in sync manner.

        :param query_ids: List of query ids (original cell ids from the input file).
        :param knn_response: MatchResult returned by space vector search service.

        :return: The fully-qualified name of the temporary table.
        """
        my_uuid = str(uuid.uuid4())[:8]
        temp_table_fqn = f"{settings.API_REQUEST_TEMP_TABLE_DATASET}.api_request_matches_{my_uuid}"
        table = bigquery.Table(temp_table_fqn, schema=bigquery_schemas.MATCH_CELL_RESULTS_SCHEMA)
        table.expires = datetime.now().astimezone() + timedelta(
            minutes=settings.API_REQUEST_TEMP_TABLE_DATASET_EXPIRATION
        )

        self.block_coo_matrix_db_client.create_table(table)

        rows_to_insert = []
        for i in range(0, len(knn_response.matches)):
            query_id = query_ids[i]
            for neighbor in knn_response.matches[i].neighbors:
                rows_to_insert.append(
                    {
                        "query_id": str(query_id),
                        "match_cas_cell_index": int(neighbor.cas_cell_index),
                        "match_score": float(neighbor.distance),
                    }
                )

        job_config = bigquery.LoadJobConfig(schema=bigquery_schemas.MATCH_CELL_RESULTS_SCHEMA)

        job = self.block_coo_matrix_db_client.load_table_from_json(
            json_rows=rows_to_insert, destination=temp_table_fqn, job_config=job_config
        )

        job.result()  # Wait for the job to complete

        return temp_table_fqn

    def get_neighborhood_distance_summary(
        self, cas_model: models.CASModel, match_temp_table_fqn: str
    ) -> t.List[schemas.QueryCellNeighborhoodCellTypeSummaryStatistics]:
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

        query_job = self.block_coo_matrix_db_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job)

    def get_neighborhood_distance_summary_dev_details(
        self, cas_model: models.CASModel, match_temp_table_fqn: str
    ) -> t.List[schemas.QueryCellNeighborhoodCellTypeSummaryStatisticsExtended]:
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

        query_job = self.block_coo_matrix_db_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job, include_extended_output=True)

    def _get_cell_metadata_by_ids(
        self, cell_ids: t.List[int], feature_names: t.List[str]
    ) -> t.List[schemas.CellariumCellMetadata]:
        """
        Get cells by CAS IDs, maintaining the order of `cell_ids`. Preserves the order of `cell_ids` in the query
        results.

        :param cell_ids: CAS cell IDs to retrieve metadata for.

        :return: List of dictionaries representing the query results, ordered according to cell_ids.
        """
        if not cell_ids:
            raise exceptions.CellMetadataDatabaseError("No cell IDs provided")

        if len(cell_ids) > settings.MAX_CELL_IDS_PER_QUERY:
            raise exceptions.CellMetadataDatabaseError(
                f"Number of cell IDs exceeds the maximum of {settings.MAX_CELL_IDS_PER_QUERY}"
            )

        with self.system_data_db_session_maker.begin() as session:
            # Creating temp table, it will be dropped automatically when the session is closed
            metadata = sa.MetaData()
            cell_info_tmp_table = sa.Table(
                f"cells_cellinfotemptable_{str(uuid.uuid4())[:8]}",
                metadata,
                sa.Column("cas_cell_index", sa.Integer),
                prefixes=["temporary"],
            )
            metadata.create_all(session.bind)

            # Insert values into the temp table
            insert_query = cell_info_tmp_table.insert().values([{"cas_cell_index": cell_id} for cell_id in cell_ids])
            session.execute(insert_query)

            # Select values from the temp table, joining with the CellInfo table
            select_columns = [getattr(models.CellInfo, feature_name) for feature_name in feature_names]
            statement = sa.select(select_columns).join(
                cell_info_tmp_table, models.CellInfo.cas_cell_index == cell_info_tmp_table.c.cas_cell_index
            )

            query_result = session.execute(statement).fetchall()

        # Reconstruct the ordered list based on cell_ids
        index_to_metadata = {row["cas_cell_index"]: row for row in query_result}
        ordered_metadata = (index_to_metadata[cell_id] for cell_id in cell_ids if cell_id in index_to_metadata)

        return [schemas.CellariumCellMetadata(**cell_metadata_row) for cell_metadata_row in ordered_metadata]

    def get_cell_metadata_by_ids(
        self, cell_ids: t.List[int], metadata_feature_names: t.List[str]
    ) -> t.List[schemas.CellariumCellMetadata]:
        """
        Get cells by CAS IDs, maintaining the order of `cell_ids`. If the number of querying cell IDs exceeds the
        maximum allowed, split the query into multiple queries. Database gets overwhelmed when querying a large number
        of cell.

        :param cell_ids: CAS cell IDs to retrieve metadata for.
        :param metadata_feature_names: List of feature names to retrieve.

        :return: List of dictionaries representing the query results, ordered according to `cell_ids`.
        """
        if len(cell_ids) >= settings.MAX_CELL_IDS_PER_QUERY:
            neighbors_metadata = []
            for i in range(0, len(cell_ids), settings.MAX_CELL_IDS_PER_QUERY):
                upper_bound_slice_index = min(len(cell_ids), i + settings.MAX_CELL_IDS_PER_QUERY)
                neighbors_metadata += self._get_cell_metadata_by_ids(
                    cell_ids=cell_ids[i:upper_bound_slice_index],
                    feature_names=metadata_feature_names,
                )
        else:
            neighbors_metadata = self._get_cell_metadata_by_ids(cell_ids=cell_ids, feature_names=metadata_feature_names)

        return neighbors_metadata
