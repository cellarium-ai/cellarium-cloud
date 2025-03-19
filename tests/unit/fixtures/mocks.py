import random
import typing as t
from unittest.mock import AsyncMock, MagicMock

import anndata
import numpy as np

from casp.services.api import data_manager, schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.db import models
from casp.services.model_inference.services import ModelInferenceService
from tests.unit.fixtures import constants


class MockCellQuotaDataManager(data_manager.CellQuotaDataManager):
    """
    Mocked version of `CellQuotaDataManager` to simulate behavior for unit tests. This mock is required because
    CellQuotaDataManager executes a SQL query with `array_agg` which is not available in sqlite.

    This mock tracks calls to its methods via `MagicMock` for verification in tests.
    """

    SQL_GET_CELLS_PROCESSED_THIS_WEEK = "SELECT mock_query"

    def __init__(self):
        """
        Initialize the mock data manager.

        The base class is prevented from performing file I/O or database setup. Instead, key behaviors
        are mocked to ensure predictable behavior in unit tests.
        """
        super().__init__()  # Safely call the base class

        self.get_cells_processed_this_week_for_user = MagicMock(side_effect=self._mock_get_cells_processed)

    def _mock_get_cells_processed(self, user: models.User) -> int:
        """
        Side effect for `get_cells_processed_this_week_for_user`.

        :param user: User instance.
        :return: Simulated number of cells processed this week.
        """
        # Simulated behavior
        if user.username == constants.USER_EMAIL_WITH_QUOTA:
            return constants.TOTAL_CELLS_PROCESSED_THIS_WEEK_FOR_USER_WITH_QUOTA

        return constants.TOTAL_CELLS_PROCESSED_THIS_WEEK_FOR_USER_WITHOUT_QUOTA


class MockMatchingClient(AsyncMock):
    """
    Class-based mock for the MatchingClient.

    This mock ensures unique neighbors are returned for each query cell during testing. The `match` method is mocked
    to simulate the behavior of retrieving nearest neighbors for a list of query embeddings. A centralized cell
    resource (`cell_source`) is used to provide valid neighbors.
    """

    def __init__(self, cell_info_data: t.List[schemas.CellariumCellMetadata], *args, **kwargs) -> None:
        """
        Initialize the MockMatchingClient with a centralized cell resource.

        The `match` method is configured as an asynchronous mock with `_match_side_effect` as its side effect.

        :param cell_info_data: A dictionary representing the centralized cell resource where keys are cell indices and
            values are associated metadata.
        :param kwargs: Additional keyword arguments passed to the parent AsyncMock.
        :param args: Additional positional arguments passed to the parent AsyncMock.
        """
        super().__init__(*args, **kwargs)
        self.cell_info_data = cell_info_data
        self.match = AsyncMock(side_effect=self._match_side_effect)

    async def _match_side_effect(self, queries: t.List[t.Any]) -> MatchResult:
        """
        Simulates the `match` method of the MatchingClient.

        For each query, it returns up to 3 unique neighbors randomly sampled from the centralized `cell_source`.

        :param queries: A list of query embeddings or identifiers for which neighbors are retrieved.

        :return: A MatchResult object containing nearest neighbors for each query.
        """
        matches = []

        # Create a mapping from cas_cell_index to CellariumCellMetadata for quick lookup
        cell_index_to_metadata = {cell.cas_cell_index: cell for cell in self.cell_info_data}

        # List of available cell indices
        cell_indices = list(cell_index_to_metadata.keys())

        for query in queries:
            num_neighbors = min(3, len(cell_indices))
            neighbors = random.sample(cell_indices, num_neighbors)

            matches.append(
                MatchResult.NearestNeighbors(
                    neighbors=[
                        MatchResult.Neighbor(
                            cas_cell_index=str(neighbor_id),
                            distance=np.random.rand(),
                        )
                        for neighbor_id in neighbors
                    ]
                )
            )
        return MatchResult(matches=matches)


class MockModelService(ModelInferenceService):
    """
    Mocked version of `ModelInferenceService` to simulate embedding behavior for unit tests.

    This mock tracks method calls and generates synthetic embeddings for input data.
    """

    def __init__(self):
        """
        Initialize the mock model service.

        The `embed_adata` method is mocked with a realistic side effect to simulate embedding generation.
        """
        super().__init__()  # Safely call the base class initializer if needed.

        # Trackable mock methods with side effects
        self.embed_adata = MagicMock(side_effect=self._mock_embed_adata)

    def _mock_embed_adata(self, adata: anndata.AnnData, model: models.CASModel) -> t.Tuple[t.List[str], np.ndarray]:
        """
        Side effect method for `embed_adata` to simulate embedding generation.

        :param adata: Anndata object containing cell data.
        :param model: CASModel instance specifying the model to use.

        :return: A tuple containing:
                 - List of observation IDs from `adata.obs_names`.
                 - Randomly generated embeddings as a NumPy array of shape (n_cells, 32).
        """
        obs_ids = adata.obs_names.tolist()
        random_embeddings = np.random.rand(len(obs_ids), 32)
        return obs_ids, random_embeddings
