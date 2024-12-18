import typing as t

import sqlalchemy as sa

from casp.data_manager import BaseDataManager
from casp.services import settings
from casp.services.api import schemas
from casp.services.api.data_manager import exceptions
from casp.services.db import models


class CellOperationsDataManager(BaseDataManager):
    """
    Data Manager for making data operations in Cellarium Cloud storage.
    """

    @classmethod
    def create_temporary_table_with_cell_ids(cls, cell_ids: t.List[int], connection: sa.engine.Connection) -> sa.Table:
        """
        Create a temporary table in postgres database and populate it with cell indexes for further joining.
        This tables acts as intermediate cache table to make receiving cell metadata query lighter and faster.

        :param cell_ids: List of cas_cell_index to insert in a temporary table
        :param connection: Database connection

        :return: An instance of :class:`sqlalchemy.Table` object
        """
        metadata = sa.MetaData(bind=connection)
        cell_info_tmp_table = sa.Table(
            "cells_cellinfotemptable",
            metadata,
            sa.Column("cas_cell_index", sa.Integer),
            prefixes=["temporary"],
            postgresql_on_commit="drop",
        )
        cell_info_tmp_table.create(bind=connection)

        # Insert values into the temp table
        rows_to_insert = [{"cas_cell_index": cell_id} for cell_id in cell_ids]
        cls.batch_bulk_insert(
            table=cell_info_tmp_table,
            rows_to_insert=rows_to_insert,
            connection=connection,
            batch_size=settings.MAX_CELL_IDS_PER_QUERY,
        )
        return cell_info_tmp_table

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
            cells_metadata = []
            for i in range(0, len(cell_ids), settings.MAX_CELL_IDS_PER_QUERY):
                upper_bound_slice_index = min(len(cell_ids), i + settings.MAX_CELL_IDS_PER_QUERY)
                cells_metadata += self._get_cell_metadata_by_ids(
                    cell_ids=cell_ids[i:upper_bound_slice_index],
                    feature_names=metadata_feature_names,
                )
        else:
            cells_metadata = self._get_cell_metadata_by_ids(cell_ids=cell_ids, feature_names=metadata_feature_names)

        return cells_metadata
