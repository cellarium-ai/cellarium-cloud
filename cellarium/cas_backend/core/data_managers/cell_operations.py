import tiledbsoma

from cellarium.cas_backend.apps.compute import schemas
from cellarium.cas_backend.core.data_managers import exceptions


class CellOperationsDataManager:
    """
    Data Manager for retrieving cell metadata from a TileDB SOMA DataFrame stored in GCS.

    The SOMA DataFrame contains the same cell metadata previously held in the Postgres ``cells_cellinfo``
    table. Using a per-model GCS-backed SOMA store removes the shared Postgres bottleneck and lets Cloud
    Run scale horizontally without a centralised DB connection limit.
    """

    def get_cell_metadata_by_ids(
        self,
        cell_metadata_uri: str,
        cell_ids: list[int],
        metadata_feature_names: list[str],
    ) -> list[schemas.CellariumCellMetadata]:
        """
        Get cell metadata by CAS cell IDs from the SOMA DataFrame, maintaining the order of ``cell_ids``.

        :param cell_metadata_uri: GCS URI pointing to the TileDB SOMA DataFrame for this model's cell metadata.
        :param cell_ids: CAS cell IDs to retrieve metadata for.
        :param metadata_feature_names: Column names to retrieve from the SOMA DataFrame.  Must be a subset
            of the fields defined on :class:`~cellarium.cas_backend.apps.compute.schemas.CellariumCellMetadata`.
            ``cas_cell_index`` is always included regardless of whether it appears in this list.

        :return: List of :class:`~cellarium.cas_backend.apps.compute.schemas.CellariumCellMetadata` objects
            ordered according to ``cell_ids``.

        :raises exceptions.CellMetadataDatabaseError: If ``cell_ids`` is empty.
        """
        if not cell_ids:
            raise exceptions.CellMetadataDatabaseError("No cell IDs provided")

        # cas_cell_index is stored as soma_joinid (the SOMA index column).
        # soma_joinid must be explicitly included in column_names — tiledbsoma does not
        # return it automatically when column_names is restricted.
        columns_to_read = [col for col in metadata_feature_names if col != "cas_cell_index"]
        soma_columns_to_read = ["soma_joinid"] + columns_to_read

        with tiledbsoma.open(cell_metadata_uri) as soma_df:
            result_table = soma_df.read(
                coords=(cell_ids,),
                column_names=soma_columns_to_read,
            ).concat()

        result_dict = result_table.to_pydict()
        # soma_joinid holds the cas_cell_index values; rename it for schema compatibility.
        index_to_row = {
            soma_joinid: {"cas_cell_index": soma_joinid, **{col: result_dict[col][i] for col in columns_to_read}}
            for i, soma_joinid in enumerate(result_dict["soma_joinid"])
        }

        return [
            schemas.CellariumCellMetadata(**index_to_row[cell_id]) for cell_id in cell_ids if cell_id in index_to_row
        ]
