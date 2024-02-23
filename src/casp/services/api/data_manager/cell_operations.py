import typing as t
import uuid
from datetime import datetime, timedelta

from google.cloud import bigquery

from casp.data_manager import BaseDataManager, sql
from casp.services import settings
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.data_manager import bigquery_response_parsers, bigquery_schemas, cellarium_general
from casp.services.db.models import CASModel, CellInfo


class CellOperationsDataManager(BaseDataManager):
    """
    Data Manager for making data operations in Cellarium Cloud storage.
    """

    # Directories for SQL templates
    CELL_ANALYSIS_TEMPLATE_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_templates/cell_operations"

    # SQL template file paths
    SQL_MATCH_METADATA = f"{CELL_ANALYSIS_TEMPLATE_DIR}/get_neighborhood_distance_summary.sql.mako"
    SQL_MATCH_METADATA_DEV_DETAILS = (
        f"{CELL_ANALYSIS_TEMPLATE_DIR}/get_neighborhood_distance_summary_dev_details.sql.mako"
    )
    SQL_GET_CELLS_BY_IDS = f"{CELL_ANALYSIS_TEMPLATE_DIR}/get_cell_metadata_by_ids.sql.mako"
    redis_cache_key_prefix = "cell_operations"

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
        self, cas_model: CASModel, match_temp_table_fqn: str
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

        query_job = self.block_coo_matrix_db_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job)

    def get_neighborhood_distance_summary_dev_details(
        self, cas_model: CASModel, match_temp_table_fqn: str
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

        query_job = self.block_coo_matrix_db_client.query(query=sql_query)

        return bigquery_response_parsers.parse_match_query_job(query_job=query_job, include_dev_details=True)

    def get_cell_metadata_by_ids(self, cell_ids: t.List[int]) -> t.List[t.Dict[str, t.Any]]:
        """
        Get cells by ids, maintaining the order of cell_ids. Try to get cell metadata from cache first, then from the
        database if not found in cache.

        :param cell_ids: Cas cell indexes to retrieve metadata for.

        :return: List of dictionaries representing the query results, ordered according to cell_ids.
        """
        # Attempt to get cell metadata from cache
        cache_keys = [f"cell_metadata:{cell_id}" for cell_id in cell_ids]
        cached_metadata_dicts = self.cache.get_many(cache_keys)
        found_ids = [int(key.split(":")[-1]) for key in cached_metadata_dicts.keys()]
        missed_ids = list(set(cell_ids) - set(found_ids))

        # Initialize a dict to hold the final result in order
        ordered_metadata = {cell_id: None for cell_id in cell_ids}

        # Update the ordered dict with cached results
        for cell_id in found_ids:
            ordered_metadata[cell_id] = cached_metadata_dicts[f"cell_metadata:{cell_id}"]

        # If there are missed_ids, retrieve them from the database
        if missed_ids:
            columns = [column.key for column in CellInfo.__table__.columns if column.key != "obs_metadata_extra"]
            with_entity_items = [getattr(CellInfo, name) for name in columns]
            database_cells = (
                self.system_data_db_session.query(CellInfo)
                .with_entities(*with_entity_items)
                .filter(CellInfo.cas_cell_index.in_(missed_ids))
                .all()
            )

            # Convert database result to dictionary format
            database_cells_dicts = [{columns[i]: value for i, value in enumerate(row)} for row in database_cells]

            # Cache newly retrieved items and update the ordered dict
            new_cache_items = {}
            for cell_dict in database_cells_dicts:
                cell_id = cell_dict["cas_cell_index"]
                cache_key = f"cell_metadata:{cell_id}"
                new_cache_items[cache_key] = cell_dict
                ordered_metadata[cell_id] = cell_dict

            # Cache the new items with a timeout
            self.cache.set_many(new_cache_items, timeout=240)

        # Extract the ordered results, filtering out any None values in case of missing IDs
        return [ordered_metadata[cell_id] for cell_id in cell_ids if ordered_metadata[cell_id] is not None]
