import json
import typing as t

import pytest

from casp.services.api import services
from casp.services.api.services.cell_operations_service import CellOperationsService
from casp.services.api.services.annotation_engines import CellOntologyResource
from casp.services.model_inference.services import RepresentationModelInferenceService
from tests.unit.fixtures import mocks


@pytest.fixture
def mock_cell_quota_service() -> services.CellQuotaService:
    """
    Pytest fixture to provide a mocked `CellQuotaDataManager`.

    :return: An instance of `MockCellQuotaDataManager`.
    """
    return services.CellQuotaService(cell_quota_dm=mocks.MockCellQuotaDataManager())


@pytest.fixture
def mock_model_service() -> mocks.MockRepresentationModelService:
    """
    Fixture for a mocked `ModelInferenceService` with embedding behavior.

    :return: A mocked `ModelInferenceService` instance that generates random embeddings.
    """
    return mocks.MockRepresentationModelService()


@pytest.fixture
def cell_operations_service_with_mocks(
    mock_cell_quota_service: services.CellQuotaService,
    mock_model_service: RepresentationModelInferenceService,
) -> CellOperationsService:
    """
    Fixture to provide a `CellOperationsService` with all mocked dependencies.

    :param mock_cell_quota_service: Mocked `CellQuotaService`.
    :param mock_model_service: Mocked `ModelInferenceService`.

    :return: An instance of `CellOperationsService` with all dependencies mocked.
    """
    service = CellOperationsService(
        cell_quota_service=mock_cell_quota_service,
        representation_model_service=mock_model_service,
    )

    return service


def load_ontology_resource_from_file() -> t.Dict[str, t.Any]:
    """
    Load the ontology resource from a predefined JSON file.

    :return: A dictionary representing the simplified ontology resource.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


@pytest.fixture
def cell_ontology_resource_mock() -> CellOntologyResource:
    """
    Fixture to provide a mocked `CellOntologyResource` instance.

    :return: A `CellOntologyResource` mock instance.
    """
    resource_dict = load_ontology_resource_from_file()
    return CellOntologyResource(cell_ontology_resource_dict=resource_dict)
