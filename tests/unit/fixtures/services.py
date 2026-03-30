import json
import typing as t
from unittest.mock import patch

import pytest

from cellarium.cas_backend.apps.compute import schemas, services
from cellarium.cas_backend.core.data_managers.cell_operations import CellOperationsDataManager
from tests.unit.fixtures import mocks

if t.TYPE_CHECKING:
    from cellarium.cas_backend.apps.compute.services.cell_operations_service import CellOperationsService
    from cellarium.cas_backend.apps.compute.services.consensus_engine.strategies.ontology_aware import (
        CellOntologyResource,
    )


@pytest.fixture
def mock_cell_quota_service() -> services.CellQuotaService:
    """
    Pytest fixture to provide a mocked `CellQuotaDataManager`.

    :return: An instance of `MockCellQuotaDataManager`.
    """
    return services.CellQuotaService(cell_quota_dm=mocks.MockCellQuotaDataManager())


@pytest.fixture
def mock_model_service() -> mocks.MockModelService:
    """
    Fixture for a mocked `ModelInferenceService` with embedding behavior.

    :return: A mocked `ModelInferenceService` instance that generates random embeddings.
    """
    return mocks.MockModelService()


@pytest.fixture
def cell_operations_service_with_mocks(
    mock_cell_quota_service: services.CellQuotaService,
    mock_model_service: mocks.MockModelService,
) -> "CellOperationsService":
    """
    Fixture to provide a `CellOperationsService` with all mocked dependencies.

    :param mock_cell_quota_service: Mocked `CellQuotaService`.
    :param mock_model_service: Mocked `ModelInferenceService`.

    :return: An instance of `CellOperationsService` with all dependencies mocked.
    """
    from cellarium.cas_backend.apps.compute.services.cell_operations_service import CellOperationsService

    service = CellOperationsService(
        cell_quota_service=mock_cell_quota_service,
        model_service=mock_model_service,
    )

    return service


@pytest.fixture
def patch_cell_operations_dm(
    cell_info_data: list[schemas.CellariumCellMetadata],
) -> t.Generator[None, None, None]:
    """
    Patch CellOperationsDataManager.get_cell_metadata_by_ids to return in-memory cell_info_data,
    avoiding real TileDB/GCS reads in unit tests.
    """

    def _mock_get(self, cell_ids: list[int], metadata_feature_names: list[str]):
        id_set = set(cell_ids)
        return [c for c in cell_info_data if c.cas_cell_index in id_set]

    with patch.object(CellOperationsDataManager, "get_cell_metadata_by_ids", _mock_get):
        yield


def load_ontology_resource_from_file() -> dict[str, t.Any]:
    """
    Load the ontology resource from a predefined JSON file.

    :return: A dictionary representing the simplified ontology resource.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


@pytest.fixture
def cell_ontology_resource_mock() -> "CellOntologyResource":
    """
    Fixture to provide a mocked `CellOntologyResource` instance.

    :return: A `CellOntologyResource` mock instance.
    """
    from cellarium.cas_backend.apps.compute.services.consensus_engine.strategies.ontology_aware import (
        CellOntologyResource,
    )

    resource_dict = load_ontology_resource_from_file()
    return CellOntologyResource(cell_ontology_resource_dict=resource_dict)
