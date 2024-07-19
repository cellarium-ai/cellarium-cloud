"""
This test module verifies the functionality of the ConsensusEngine ontology-aware method. It includes loading
ontology resources, processing match results, and comparing the output against expected responses.
`cell_type_count` method isn't covered as current version is deprecated and has to be rewritten with use of PostgreSQL
instead of BigQuery as it is now.

The tests make use of fixtures to provide necessary mock data and configurations.
"""

import json
import typing as t
from unittest.mock import Mock

import pytest

from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.services import consensus_engine


def load_ontology_resource_from_file() -> t.Dict[str, t.Any]:
    """
    Load the ontology resource from a predefined JSON file.

    :return: A dictionary representing the simplified ontology resource.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


def load_expected_response() -> list[schemas.QueryCellNeighborhoodOntologyAware]:
    """
    Load the expected response from a predefined JSON file.

    :return: A list of `QueryCellAnnotationOntologyAware` instances.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/ontology_aware_expected_response.json"

    with open(filepath, "r") as file:
        json_data = json.load(file)

    return [schemas.QueryCellNeighborhoodOntologyAware(**item) for item in json_data]


@pytest.fixture
def cell_ontology_resource_mock() -> consensus_engine.CellOntologyResource:
    """
    Fixture to provide a mocked `CellOntologyResource` instance.

    :return: A `CellOntologyResource` mock instance.
    """
    resource_dict = load_ontology_resource_from_file()
    return consensus_engine.CellOntologyResource(cell_ontology_resource_dict=resource_dict)


@pytest.fixture
def cell_operations_dm_mock() -> Mock:
    """
    Fixture to provide a mocked `CellOperationsDataManager` instance.

    :return: A `CellOperationsDataManager` mock instance.
    """
    dm_mock = Mock()
    dm_mock.get_cell_metadata_by_ids.return_value = [
        schemas.CellariumCellMetadata(
            cas_cell_index=48033450,
            cell_type="alpha-beta T cell",
            cell_type_ontology_term_id="CL:0000789",
        ),
        schemas.CellariumCellMetadata(
            cas_cell_index=1384010152, cell_type="T cell", cell_type_ontology_term_id="CL:0000084"
        ),
    ]
    return dm_mock


@pytest.fixture
def consensus_engine_mock(cell_ontology_resource_mock, cell_operations_dm_mock) -> consensus_engine.ConsensusEngine:
    """
    Fixture to provide a `ConsensusEngine` instance with mocked dependencies.

    :param cell_ontology_resource_mock: A fixture providing a mocked `CellOntologyResource` instance.
    :param cell_operations_dm_mock: A fixture providing a mocked `CellOperationsDataManager` instance.
    :return:
    """
    strategy = consensus_engine.CellTypeOntologyAwareConsensusStrategy(
        prune_threshold=0.1,
        weighting_prefactor=1.0,
        cell_ontology_resource=cell_ontology_resource_mock,
        cell_operations_dm=cell_operations_dm_mock,
    )
    return consensus_engine.ConsensusEngine(strategy=strategy)


@pytest.fixture
def knn_query_result_mock() -> MatchResult:
    """
    Mocked KNN query result data to simulate neighbor matching results.

    :return: A `MatchResult` instance with mocked neighbor matching results.
    """
    return MatchResult(
        matches=[
            MatchResult.NearestNeighbors(
                neighbors=[
                    MatchResult.Neighbor(cas_cell_index="48033450", distance=0.9),
                    MatchResult.Neighbor(cas_cell_index="1384010152", distance=0.87),
                ]
            ),
        ]
    )


def test_summarize_query_neighbor_context_ontology_aware(consensus_engine_mock, knn_query_result_mock):
    """
    Test the ontology-aware query neighbor context summarization.

    Verifies that the `summarize_query_neighbor_context` method of the `ConsensusEngine`
    returns correct summarizations when using the ontology-aware method. It checks
    that the output matches an expected response loaded from a predefined JSON file, ensuring
    that the summarization logic correctly processes and weights the match results based on
    cell type ontologies.

    Assertions:
    - Asserts that the output is a list of `QueryCellAnnotationOntologyAware` instances.
    - Asserts that the output matches the expected response defined in `expected_response.json`.

    :param consensus_engine_mock: A `ConsensusEngine` instance with mocked dependencies.
    :param knn_query_result_mock: Mocked `MatchResult` data to simulate neighbor matching results.
    """
    query_ids = ["query_ACGTGCTGATG_1"]

    result = consensus_engine_mock.summarize(query_ids=query_ids, knn_query=knn_query_result_mock)
    expected = load_expected_response()
    # Assertions to validate the output
    assert isinstance(result, list), "The result should be a list"
    assert all(
        isinstance(item, schemas.QueryCellNeighborhoodOntologyAware) for item in result
    ), "All items should be of type QueryCellAnnotationOntologyAware"

    assert expected == result, "The expected response should match the result"
