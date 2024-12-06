import typing as t
from datetime import datetime

import sqlalchemy as sa

from casp.data_manager import BaseDataManager
from casp.services import settings
from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.data_manager import exceptions
from casp.services.db import models


class CellOperationsDataManager(BaseDataManager):
    """
    Data Manager for making data operations in Cellarium Cloud storage.
    """

    @staticmethod
    def create_temporary_table_with_cell_ids(cell_ids: t.List[int], connection: sa.engine.Connection) -> sa.Table:
        metadata = sa.MetaData(bind=connection)
        cell_info_tmp_table = sa.Table(
            f"cells_cellinfotemptable",
            metadata,
            sa.Column("cas_cell_index", sa.Integer),
            prefixes=["temporary"],
        )
        cell_info_tmp_table.create(bind=connection)

        # Insert values into the temp table
        insert_query = cell_info_tmp_table.insert().values([{"cas_cell_index": cell_id} for cell_id in cell_ids])
        connection.execute(insert_query)
        return cell_info_tmp_table

    @classmethod
    def create_temporary_table_with_knn_match_data(
        cls, query_ids: t.List[str], knn_response: MatchResult, connection: sa.engine.Connection
    ) -> sa.Table:
        metadata = sa.MetaData(bind=connection)
        # metadata = sa.MetaData(bind=session.bind)
        match_tmp_table = sa.Table(
            f"cells_requestmatchtemptable",
            metadata,
            sa.Column("query_id", sa.String, nullable=False),
            sa.Column("match_cas_cell_index", sa.Integer, nullable=False),
            sa.Column("match_score", sa.Float, nullable=False),
            sa.Column("created_at", sa.TIMESTAMP, nullable=False, default=datetime.utcnow),
            prefixes=["temporary"],
        )
        match_tmp_table.create(bind=connection)

        rows_to_insert = []
        for query_id, match in zip(query_ids, knn_response.matches):
            for neighbor in match.neighbors:
                rows_to_insert.append(
                    {
                        "query_id": str(query_id),
                        "match_cas_cell_index": int(neighbor.cas_cell_index),
                        "match_score": float(neighbor.distance),
                        "created_at": datetime.utcnow(),
                    }
                )

        # Insert data in batches because this query becomes large
        cls.batch_bulk_insert(
            table=match_tmp_table, rows_to_insert=rows_to_insert, connection=connection, batch_size=20000
        )

        return match_tmp_table

    def _get_cell_metadata_by_ids(
        self,
        cell_ids: t.List[int],
        feature_names: t.List[str],
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
            connection = session.connection()
            cell_info_tmp_table = self.create_temporary_table_with_cell_ids(cell_ids=cell_ids, connection=connection)

            # Select values from the temp table, joining with the CellInfo table
            select_columns = [getattr(models.CellInfo, feature_name) for feature_name in feature_names]
            statement = sa.select(select_columns).join(
                cell_info_tmp_table, models.CellInfo.cas_cell_index == cell_info_tmp_table.c.cas_cell_index
            )

            query_result = connection.execute(statement).fetchall()

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

    @staticmethod
    def parse_query_statistics_from_db_results(
        db_results: t.Iterable[sa.engine.Row],
    ) -> t.List[schemas.QueryCellNeighborhoodCellTypeSummaryStatistics]:
        query_cell_statistics = []

        for row in db_results:
            # Each row contains `query_id` and `matches` (list of JSON objects)
            query_cell_id = row["query_cell_id"]
            matches = [
                schemas.NeighborhoodCellTypeSummaryStatistics(**match)
                for match in row["matches"]  # Assuming `matches` is already a JSON array
            ]

            # Create the Pydantic model for each query
            query_cell_statistics.append(
                schemas.QueryCellNeighborhoodCellTypeSummaryStatistics(
                    query_cell_id=query_cell_id,
                    matches=matches,
                )
            )

        return query_cell_statistics

    @staticmethod
    def calculate_query_cell_neighborhood_summary_statistics(
        knn_match_data_table: sa.Table,
        knn_unique_neighbors_cell_info_index_table: sa.Table,
        connection: sa.engine.Connection,
    ) -> t.List[schemas.QueryCellNeighborhoodCellTypeSummaryStatistics]:
        # Define the subquery to compute and order the aggregate statistics
        ordered_subquery = (
            sa.select(
                knn_match_data_table.c.query_id.label("query_cell_id"),  # Rename query_id to query_cell_id
                models.CellInfo.cell_type,
                sa.func.count().label("cell_count"),
                sa.func.min(knn_match_data_table.c.match_score).label("min_distance"),
                sa.func.max(knn_match_data_table.c.match_score).label("max_distance"),
                sa.func.percentile_cont(0.25).within_group(knn_match_data_table.c.match_score).label("p25_distance"),
                sa.func.percentile_cont(0.5).within_group(knn_match_data_table.c.match_score).label("median_distance"),
                sa.func.percentile_cont(0.75).within_group(knn_match_data_table.c.match_score).label("p75_distance"),
            )
            .select_from(
                knn_match_data_table.join(
                    knn_unique_neighbors_cell_info_index_table,
                    knn_match_data_table.c.match_cas_cell_index
                    == knn_unique_neighbors_cell_info_index_table.c.cas_cell_index,
                ).join(
                    models.CellInfo,
                    knn_match_data_table.c.match_cas_cell_index == models.CellInfo.cas_cell_index,
                )
            )
            .group_by(
                knn_match_data_table.c.query_id,
                models.CellInfo.cell_type,
            )
            .order_by(
                knn_match_data_table.c.query_id,  # Group by query_id and order results
                sa.desc(sa.func.count()),  # Ensure rows are ordered by cell_count descending
            )
            .subquery()
        )
        # Use the ordered subquery to aggregate the results into JSON objects
        query = (
            sa.select(
                ordered_subquery.c.query_cell_id,  # Use the renamed column
                sa.func.array_agg(
                    sa.func.json_build_object(
                        "cell_type",
                        ordered_subquery.c.cell_type,
                        "cell_count",
                        ordered_subquery.c.cell_count,
                        "min_distance",
                        ordered_subquery.c.min_distance,
                        "max_distance",
                        ordered_subquery.c.max_distance,
                        "p25_distance",
                        ordered_subquery.c.p25_distance,
                        "median_distance",
                        ordered_subquery.c.median_distance,
                        "p75_distance",
                        ordered_subquery.c.p75_distance,
                    )
                ).label("matches"),
            )
            .group_by(ordered_subquery.c.query_cell_id)  # Group by the renamed column
            .order_by(ordered_subquery.c.query_cell_id)  # Order by the renamed column
        )

        # Execute the query using the connection
        results = connection.execute(query).fetchall()

        return CellOperationsDataManager.parse_query_statistics_from_db_results(db_results=results)
