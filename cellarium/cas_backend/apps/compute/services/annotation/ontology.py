"""
Shared cell-ontology scoring primitives used by both the kNN consensus path
(CellTypeOntologyAwareConsensusStrategy) and the classification uplift path
(ClassificationOntologyUpliftStrategy). The propagation maths must not be duplicated.
"""

import json
import typing as t

from smart_open import open

from cellarium.cas_backend.apps.compute import schemas


class CellOntologyResource:
    """
    Handles cell ontology resource data, such as the ancestors dictionary and cell ontology term ID mappings.

    Attributes:
        ancestors_dictionary: Maps cell ontology term IDs to their ancestor IDs.
        ontology_term_id_to_name_dict: Maps cell ontology term IDs to cell type names.

    Raises:
        ValueError: If either ``ancestors_dictionary`` or ``cell_ontology_term_id_to_cell_type`` is missing from the
            provided dictionary.
    """

    def __init__(
        self,
        ontology_resource_uri: str | None = None,
        cell_ontology_resource_dict: dict[str, t.Any] | None = None,
    ):
        if cell_ontology_resource_dict is None:
            if ontology_resource_uri is None:
                raise ValueError(
                    "`ontology_resource_uri` must be provided when `cell_ontology_resource_dict` is not given"
                )
            with open(ontology_resource_uri, "rb") as f:
                cell_ontology_resource_dict = json.loads(f.read())

        if "ancestors_dictionary" not in cell_ontology_resource_dict:
            raise ValueError("`ancestors_dictionary` is not found in the cell ontology resource dictionary")
        if "cell_ontology_term_id_to_cell_type" not in cell_ontology_resource_dict:
            raise ValueError(
                "`cell_ontology_term_id_to_cell_type` mapping is not found in the cell ontology resource dictionary"
            )

        self.ancestors_dictionary = cell_ontology_resource_dict["ancestors_dictionary"]
        self.ontology_term_id_to_name_dict = cell_ontology_resource_dict["cell_ontology_term_id_to_cell_type"]
        self.children_dictionary: dict[str, list[str]] | None = cell_ontology_resource_dict.get("children_dictionary")
        self.shortest_path_lengths_from_cell_root: dict[str, int] | None = cell_ontology_resource_dict.get(
            "shortest_path_lengths_from_cell_root"
        )
        self.longest_path_lengths_from_cell_root: dict[str, int] | None = cell_ontology_resource_dict.get(
            "longest_path_lengths_from_cell_root"
        )

    @property
    def cl_names(self) -> list[str]:
        return sorted(self.ontology_term_id_to_name_dict.keys())


def normalize_cl_term_id(term_id: str) -> str:
    """
    Normalize a cell ontology term id to the colon-separated spelling the ontology resource uses.

    :param term_id: Cell ontology term id, e.g. ``CL:0000789`` or ``CL_0000789``.

    :return: ``term_id`` with a leading ``CL_`` rewritten to ``CL:``; unchanged otherwise.
    """
    if term_id.startswith("CL_"):
        return term_id.replace("CL_", "CL:", 1)
    return term_id


def accumulate_ontology_scores(
    *, resource: CellOntologyResource, term_ids: t.Sequence[str], weights: t.Sequence[float]
) -> tuple[dict[str, float], int]:
    """
    Accumulate weights onto cell ontology terms and propagate them to every transitive ancestor.

    This is the shared propagation primitive for both the kNN consensus path (neighbour cell types
    weighted by distance) and the classification uplift path (model class labels weighted by softmax
    probability): the maths is identical for both callers and must not be duplicated.

    :param resource: Cell ontology resource supplying the ancestor graph.
    :param term_ids: Cell ontology term ids to accumulate, one per weight.
    :param weights: Weight to add for each term id, one per term id.

    :return: A tuple of the accumulated per-term scores (seeded at ``0.0`` for every term in
        ``resource.ancestors_dictionary``, in that dict's key order) and the count of term ids not found
        in the resource.
    """
    scores = dict.fromkeys(resource.ancestors_dictionary.keys(), 0.0)
    unrecognized_count = 0

    for term_id, weight in zip(term_ids, weights, strict=False):
        if term_id not in scores:
            unrecognized_count += 1
            continue
        scores[term_id] += weight
        for ancestor in resource.ancestors_dictionary[term_id]:
            scores[ancestor] += weight

    return scores, unrecognized_count


def build_ontology_matches(
    *, resource: CellOntologyResource, scores: dict[str, float], prune_threshold: float
) -> list[schemas.NeighborhoodSummaryOntologyAware]:
    """
    Prune low-scoring ontology terms and project the remainder into response schema objects.

    :param resource: Cell ontology resource supplying the CL id to cell-type name mapping.
    :param scores: Per-term accumulated scores, as returned by :func:`accumulate_ontology_scores`.
    :param prune_threshold: Terms scoring below this value are dropped. No filtering when ``0.0``.

    :return: One :class:`~cellarium.cas_backend.apps.compute.schemas.NeighborhoodSummaryOntologyAware` per
        surviving term, in ``scores`` iteration order.
    """
    if prune_threshold > 0.0:
        scores = {k: v for k, v in scores.items() if v >= prune_threshold}

    return [
        schemas.NeighborhoodSummaryOntologyAware(
            score=score,
            cell_type_ontology_term_id=cell_type_ontology_term_id,
            cell_type=resource.ontology_term_id_to_name_dict[cell_type_ontology_term_id],
        )
        for cell_type_ontology_term_id, score in scores.items()
    ]
