"""
This test module verifies ClassificationOntologyUpliftStrategy: uplifting a classification model's
per-cell class probabilities over the cell type ontology, and comparing the result against a golden
expected response.
"""

import json
import typing as t

import numpy as np
import pytest

from cellarium.cas_backend.apps.compute import schemas
from cellarium.cas_backend.apps.compute.services import exceptions
from cellarium.cas_backend.apps.compute.services.annotation.classification import ClassificationOntologyUpliftStrategy
from cellarium.cas_backend.apps.compute.services.annotation.ontology import (
    CellOntologyResource,
    accumulate_ontology_scores,
)


def load_ontology_resource_from_file() -> dict[str, t.Any]:
    """
    Load the ontology resource from a predefined JSON file.

    :return: A dictionary representing the simplified ontology resource.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


def load_expected_response() -> list[schemas.QueryCellNeighborhoodOntologyAware]:
    """
    Load the expected classification response from a predefined JSON file.

    :return: A list of `QueryCellNeighborhoodOntologyAware` instances.
    """
    filepath = "tests/unit/test_classification_fixtures/classification_expected_response.json"
    with open(filepath) as file:
        json_data = json.load(file)

    return [schemas.QueryCellNeighborhoodOntologyAware(**item) for item in json_data]


@pytest.fixture
def cell_ontology_resource() -> CellOntologyResource:
    """
    Fixture to provide a `CellOntologyResource` instance built from the mini ontology fixture.

    :return: A `CellOntologyResource` instance.
    """
    resource_dict = load_ontology_resource_from_file()
    return CellOntologyResource(cell_ontology_resource_dict=resource_dict)


@pytest.fixture
def strategy(cell_ontology_resource: CellOntologyResource) -> ClassificationOntologyUpliftStrategy:
    """
    Fixture to provide a `ClassificationOntologyUpliftStrategy` instance with a `0.1` prune threshold.

    :param cell_ontology_resource: A fixture providing a `CellOntologyResource` instance.
    :return: A `ClassificationOntologyUpliftStrategy` instance.
    """
    return ClassificationOntologyUpliftStrategy(cell_ontology_resource=cell_ontology_resource, prune_threshold=0.1)


def test_summarize_uplifts_probabilities_over_ontology(strategy: ClassificationOntologyUpliftStrategy):
    """
    `summarize` propagates each label's probability onto itself and every transitive ancestor, matching a
    golden expected response computed independently.

    Assertions:
    - The result matches the expected response defined in `classification_expected_response.json`.
    - Ancestor terms score at least as high as the leaf label that was actually predicted.
    - Matches come back in `ancestors_dictionary` key order, not sorted by score.
    """
    query_cell_ids = ["query_cell_1", "query_cell_2"]
    labels = ["CL:0000789", "CL:0000084", "CL:9999999"]
    probabilities = np.array(
        [
            [0.8437947344813395, 0.11419519938459449, 0.04201006613406605],
            [0.04527850074362907, 0.9094429985127419, 0.04527850074362907],
        ]
    )

    result = strategy.summarize(query_cell_ids=query_cell_ids, probabilities=probabilities, labels=labels)
    expected = load_expected_response()

    assert isinstance(result, list)
    assert all(isinstance(item, schemas.QueryCellNeighborhoodOntologyAware) for item in result)
    assert len(result) == len(expected)

    for result_item, expected_item in zip(result, expected, strict=False):
        assert result_item.query_cell_id == expected_item.query_cell_id
        assert result_item.total_neighbors == expected_item.total_neighbors
        assert result_item.total_neighbors_unrecognized == expected_item.total_neighbors_unrecognized
        assert result_item.total_weight == pytest.approx(expected_item.total_weight)
        assert len(result_item.matches) == len(expected_item.matches)
        for result_match, expected_match in zip(result_item.matches, expected_item.matches, strict=False):
            assert result_match.cell_type == expected_match.cell_type
            assert result_match.cell_type_ontology_term_id == expected_match.cell_type_ontology_term_id
            assert result_match.score == pytest.approx(expected_match.score)


def test_summarize_raises_when_all_labels_unknown(strategy: ClassificationOntologyUpliftStrategy):
    """`summarize` raises `ClassificationLabelSpaceError` when no label maps onto the ontology resource."""
    with pytest.raises(exceptions.ClassificationLabelSpaceError):
        strategy.summarize(
            query_cell_ids=["query_cell_1"],
            probabilities=np.array([[0.5, 0.5]]),
            labels=["CL:9999999", "CL:9999998"],
        )


def test_accumulate_ontology_scores_matches_socam_propagation():
    """
    Cross-validate ``accumulate_ontology_scores`` against SOCAM's ``propagate_probs`` math.

    SOCAM propagates raw softmax probabilities within its trained label set via matrix
    multiplication:

        propagated[n, k] = sum_c( descendant_tensor[k, c] * probs[n, c] )

    where ``descendant_tensor[k, c] = 1`` if c is a descendant of k, with the diagonal set
    to 1 so each node is its own descendant (self-inclusive).

    ``accumulate_ontology_scores`` expresses the same computation from the opposite direction:
    for each trained label c, add its weight to c itself and to every transitive ancestor.

    The two formulations are algebraically equivalent.  This test proves that explicitly:
    it builds a toy 4-node ontology, constructs SOCAM's descendant tensor from the same
    ancestor relationships, runs both propagations, and asserts identical scores for every
    node in the trained set across two cells with different probability distributions.
    """
    # -----------------------------------------------------------------------
    # Toy ontology: root → lymphocyte → T cell → alpha-beta T cell (linear chain)
    # -----------------------------------------------------------------------
    cl_root = "CL:0000000"
    cl_lymphocyte = "CL:0000542"
    cl_t_cell = "CL:0000084"
    cl_ab_t_cell = "CL:0000789"
    trained_labels = [cl_root, cl_lymphocyte, cl_t_cell, cl_ab_t_cell]

    # ancestors_dictionary: each term maps to its transitive ancestors (not including self),
    # which is the format CellOntologyResource expects.
    ancestors_dict = {
        cl_root: [],
        cl_lymphocyte: [cl_root],
        cl_t_cell: [cl_root, cl_lymphocyte],
        cl_ab_t_cell: [cl_root, cl_lymphocyte, cl_t_cell],
    }
    resource = CellOntologyResource(
        cell_ontology_resource_dict={
            "ancestors_dictionary": ancestors_dict,
            "cell_ontology_term_id_to_cell_type": {t: t for t in trained_labels},
        }
    )

    # -----------------------------------------------------------------------
    # Build SOCAM's descendant_tensor_cc from the same relationships.
    # descendant_tensor[k, c] = 1 if c is a descendant of k (diagonal = 1: self-inclusive).
    # -----------------------------------------------------------------------
    n = len(trained_labels)
    idx = {label: i for i, label in enumerate(trained_labels)}
    descendant_tensor = np.zeros((n, n), dtype=np.float64)
    for label in trained_labels:
        descendant_tensor[idx[label], idx[label]] = 1.0  # self
        for ancestor in ancestors_dict[label]:
            descendant_tensor[idx[ancestor], idx[label]] = 1.0  # ancestor row, descendant col

    # -----------------------------------------------------------------------
    # Two test cells with deliberately different probability distributions.
    # -----------------------------------------------------------------------
    raw_probs = np.array(
        [
            [0.05, 0.10, 0.15, 0.70],  # cell 0: mostly alpha-beta T cell
            [0.40, 0.30, 0.20, 0.10],  # cell 1: spread across all labels
        ],
        dtype=np.float64,
    )

    # SOCAM propagation: propagated[n, k] = sum_c( descendant_tensor[k, c] * probs[n, c] )
    socam_propagated = np.einsum("nc,kc->nk", raw_probs, descendant_tensor)

    # CAS propagation via accumulate_ontology_scores, one cell at a time.
    for cell_idx in range(len(raw_probs)):
        scores, unrecognized = accumulate_ontology_scores(
            resource=resource,
            term_ids=trained_labels,
            weights=raw_probs[cell_idx],
        )

        assert unrecognized == 0, "all trained labels must be recognised by the resource"

        for label in trained_labels:
            assert scores[label] == pytest.approx(socam_propagated[cell_idx, idx[label]]), (
                f"cell {cell_idx}, label {label!r}: "
                f"CAS={scores[label]:.6f}, SOCAM={socam_propagated[cell_idx, idx[label]]:.6f}"
            )
