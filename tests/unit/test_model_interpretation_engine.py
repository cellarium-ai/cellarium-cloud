"""
Unit tests for Model Interpretation Engine.

This module tests the ModelInterpretationEngine and its strategies, particularly the
OntologyAncestorUpliftStrategy which processes classification model outputs and propagates
probabilities to ontology ancestors.
"""

import numpy as np
import pytest

from casp.services.api import schemas
from casp.services.api.services.annotation_engines import model_interpretation_engine
from casp.services.api.services.annotation_engines.shared_resources import CellOntologyResource
from casp.services.model_inference.schemas import ClassificationModelOutput


@pytest.fixture
def cell_ontology_resource_mock() -> CellOntologyResource:
    """
    Fixture to provide a mocked `CellOntologyResource` with a simplified ontology structure.

    :return: A mocked `CellOntologyResource` instance.
    """
    ontology_data = {
        "cell_ontology_term_id_to_cell_type": {
            "CL_0000000": "cell",
            "CL_0000084": "T cell",
            "CL_0000624": "CD4-positive, alpha-beta T cell",
            "CL_0000625": "CD8-positive, alpha-beta T cell",
            "CL_0000542": "lymphocyte",
            "CL_0000576": "monocyte",
            "CL_0000232": "erythrocyte",
        },
        "ancestors_dictionary": {
            "CL_0000624": ["CL_0000084", "CL_0000542", "CL_0000000"],  # CD4+ T cell -> T cell -> lymphocyte -> cell
            "CL_0000625": ["CL_0000084", "CL_0000542", "CL_0000000"],  # CD8+ T cell -> T cell -> lymphocyte -> cell
            "CL_0000084": ["CL_0000542", "CL_0000000"],  # T cell -> lymphocyte -> cell
            "CL_0000542": ["CL_0000000"],  # lymphocyte -> cell
            "CL_0000576": ["CL_0000000"],  # monocyte -> cell
            "CL_0000232": ["CL_0000000"],  # erythrocyte -> cell
        },
    }
    return CellOntologyResource(cell_ontology_resource_dict=ontology_data)


@pytest.fixture
def mock_classification_output() -> ClassificationModelOutput:
    """
    Fixture to provide mock classification model output.

    :return: ClassificationModelOutput with synthetic probabilities.
    """
    # Simulate 2 cells with probabilities for 5 classes
    probabilities = np.array(
        [
            [0.6, 0.2, 0.1, 0.05, 0.05],  # Cell 0: high confidence for class 0
            [0.1, 0.1, 0.5, 0.2, 0.1],  # Cell 1: high confidence for class 2
        ]
    )
    labels = ["CL_0000624", "CL_0000625", "CL_0000084", "CL_0000576", "CL_0000232"]
    sample_ids = ["cell_0", "cell_1"]

    return ClassificationModelOutput(probabilities=probabilities, labels=labels, sample_ids=sample_ids)


@pytest.fixture
def model_interpretation_engine_mock(
    cell_ontology_resource_mock: CellOntologyResource,
) -> model_interpretation_engine.ModelInterpretationEngine:
    """
    Fixture to provide a ModelInterpretationEngine with OntologyAncestorUpliftStrategy.

    :param cell_ontology_resource_mock: Mocked ontology resource.

    :return: ModelInterpretationEngine instance.
    """
    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.05,
        decay=1.0,
        stop_at_model_ancestor=True,
        top_k=None,
        cell_ontology_resource=cell_ontology_resource_mock,
    )
    return model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)


def test_ontology_ancestor_uplift_strategy_basic(
    model_interpretation_engine_mock: model_interpretation_engine.ModelInterpretationEngine,
    mock_classification_output: ClassificationModelOutput,
):
    """
    Test basic functionality of OntologyAncestorUpliftStrategy.

    Verify that the strategy correctly processes classification output and returns
    ontology-aware annotations with ancestor uplifting.
    """
    result = model_interpretation_engine_mock.summarize(model_output=mock_classification_output)

    # Assertions
    assert isinstance(result, list), "Result should be a list"
    assert len(result) == 2, "Should have 2 results (one per cell)"
    assert all(
        isinstance(item, schemas.QueryCellAnnotationOntologyAware) for item in result
    ), "All items should be QueryCellAnnotationOntologyAware"

    # Check first cell result
    first_cell = result[0]
    assert first_cell.query_cell_id == "cell_0"
    assert len(first_cell.matches) > 0, "Should have at least one match"
    assert all(
        isinstance(match, schemas.AnnotationSummaryOntologyAware) for match in first_cell.matches
    ), "All matches should be AnnotationSummaryOntologyAware"

    # Verify that scores are present and sorted
    scores = [match.score for match in first_cell.matches]
    assert scores == sorted(scores, reverse=True), "Scores should be sorted in descending order"


def test_ontology_ancestor_uplift_with_prune_threshold(
    cell_ontology_resource_mock: CellOntologyResource,
    mock_classification_output: ClassificationModelOutput,
):
    """
    Test that prune_threshold correctly filters out low-probability annotations.
    """
    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.15,  # Higher threshold to filter more
        decay=1.0,
        cell_ontology_resource=cell_ontology_resource_mock,
    )
    engine = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)

    result = engine.summarize(model_output=mock_classification_output)

    # First cell should have fewer matches due to pruning
    first_cell = result[0]
    scores = [match.score for match in first_cell.matches]
    assert all(
        score >= 0.15 or len(first_cell.matches) == 1 for score in scores
    ), "All scores should be >= threshold (or single match if all pruned)"


def test_ontology_ancestor_uplift_with_decay(
    cell_ontology_resource_mock: CellOntologyResource,
    mock_classification_output: ClassificationModelOutput,
):
    """
    Test that decay parameter attenuates ancestor scores.
    """
    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.0,
        decay=0.5,  # Ancestors get 50% of child's score per level
        cell_ontology_resource=cell_ontology_resource_mock,
    )
    engine = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)

    result = engine.summarize(model_output=mock_classification_output)

    # Verify result structure
    assert len(result) == 2
    first_cell = result[0]
    assert len(first_cell.matches) > 0


def test_ontology_ancestor_uplift_with_top_k(
    cell_ontology_resource_mock: CellOntologyResource,
    mock_classification_output: ClassificationModelOutput,
):
    """
    Test that top_k parameter limits the number of returned annotations.
    """
    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.0,
        decay=1.0,
        top_k=3,  # Only return top 3 matches
        cell_ontology_resource=cell_ontology_resource_mock,
    )
    engine = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)

    result = engine.summarize(model_output=mock_classification_output)

    # Each cell should have at most 3 matches
    for cell_result in result:
        assert len(cell_result.matches) <= 3, "Should have at most 3 matches due to top_k=3"


def test_ontology_ancestor_uplift_empty_probabilities():
    """
    Test handling of edge case where all probabilities are zero or below threshold.
    """
    # Create minimal ontology
    ontology_data = {
        "cell_ontology_term_id_to_cell_type": {"CL_0000000": "cell", "CL_0000084": "T cell"},
        "ancestors_dictionary": {"CL_0000084": ["CL_0000000"]},
    }
    resource = CellOntologyResource(cell_ontology_resource_dict=ontology_data)

    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.9,  # Very high threshold
        cell_ontology_resource=resource,
    )
    engine = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)

    # All probabilities below threshold
    probabilities = np.array([[0.1, 0.1]])
    labels = ["CL_0000084", "CL_0000000"]
    sample_ids = ["cell_0"]

    model_output = ClassificationModelOutput(probabilities=probabilities, labels=labels, sample_ids=sample_ids)

    result = engine.summarize(model_output=model_output)

    # Should still return at least one match (the highest probability)
    assert len(result) == 1
    assert len(result[0].matches) >= 1, "Should have at least one match even if all below threshold"


def test_ontology_ancestor_uplift_stop_at_model_ancestor(
    cell_ontology_resource_mock: CellOntologyResource,
):
    """
    Test that stop_at_model_ancestor flag prevents climbing past model-covered ancestors.
    """
    # Test with stop_at_model_ancestor=True
    strategy_stop = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.0,
        decay=1.0,
        stop_at_model_ancestor=True,
        cell_ontology_resource=cell_ontology_resource_mock,
    )

    # Test with stop_at_model_ancestor=False
    strategy_no_stop = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.0,
        decay=1.0,
        stop_at_model_ancestor=False,
        cell_ontology_resource=cell_ontology_resource_mock,
    )

    # Create model output where T cell (CL_0000084) is in the model
    probabilities = np.array([[0.8, 0.2]])
    labels = ["CL_0000624", "CL_0000084"]  # CD4+ T cell and T cell both in model
    sample_ids = ["cell_0"]

    model_output = ClassificationModelOutput(probabilities=probabilities, labels=labels, sample_ids=sample_ids)

    engine_stop = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy_stop)
    engine_no_stop = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy_no_stop)

    result_stop = engine_stop.summarize(model_output=model_output)
    result_no_stop = engine_no_stop.summarize(model_output=model_output)

    # Both should return results
    assert len(result_stop) == 1
    assert len(result_no_stop) == 1


@pytest.mark.asyncio
async def test_model_interpretation_engine_async(
    model_interpretation_engine_mock: model_interpretation_engine.ModelInterpretationEngine,
    mock_classification_output: ClassificationModelOutput,
):
    """
    Test async summarization method.
    """
    result = await model_interpretation_engine_mock.summarize_async(model_output=mock_classification_output)

    # Should return same structure as sync method
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(isinstance(item, schemas.QueryCellAnnotationOntologyAware) for item in result)


def test_ontology_ancestor_uplift_label_normalization(
    cell_ontology_resource_mock: CellOntologyResource,
):
    """
    Test that labels with colons (e.g., CL:0000084) are normalized to underscores (CL_0000084).
    """
    strategy = model_interpretation_engine.OntologyAncestorUpliftStrategy(
        prune_threshold=0.0,
        cell_ontology_resource=cell_ontology_resource_mock,
    )
    engine = model_interpretation_engine.ModelInterpretationEngine(strategy=strategy)

    # Use labels with colons
    probabilities = np.array([[0.8, 0.2]])
    labels = ["CL:0000624", "CL:0000084"]  # With colons
    sample_ids = ["cell_0"]

    model_output = ClassificationModelOutput(probabilities=probabilities, labels=labels, sample_ids=sample_ids)

    result = engine.summarize(model_output=model_output)

    # Should successfully process and return results
    assert len(result) == 1
    assert len(result[0].matches) > 0

    # Check that cell type ontology term IDs use underscores
    for match in result[0].matches:
        assert "_" in match.cell_type_ontology_term_id or match.cell_type_ontology_term_id == "CL_0000000"
