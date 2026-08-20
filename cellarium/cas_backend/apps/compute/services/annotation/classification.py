"""
Classification annotation strategy: uplift a classification model's per-cell class probabilities over the
cell type ontology, producing a response shaped like the kNN consensus ontology-aware strategy.
"""

import numpy as np

from cellarium.cas_backend.apps.compute import schemas
from cellarium.cas_backend.apps.compute.services import exceptions
from cellarium.cas_backend.apps.compute.services.annotation.ontology import (
    CellOntologyResource,
    accumulate_ontology_scores,
    build_ontology_matches,
    normalize_cl_term_id,
)


class ClassificationOntologyUpliftStrategy:
    """
    Uplift a classification model's per-cell class probabilities over the cell type ontology.

    Unlike CellTypeOntologyAwareConsensusStrategy, which propagates kNN-neighbour distance weights, this
    strategy propagates a classification model's own softmax probabilities. Both paths share the
    propagation maths in :mod:`cellarium.cas_backend.apps.compute.services.annotation.ontology`.

    :param cell_ontology_resource: Cell ontology resource supplying the ancestor graph and CL id to
        cell-type name mapping.
    :param prune_threshold: Terms scoring below this value are dropped from the response. No filtering
        when ``0.0``.
    """

    def __init__(self, *, cell_ontology_resource: CellOntologyResource, prune_threshold: float):
        self.cell_ontology_resource = cell_ontology_resource
        self.prune_threshold = prune_threshold

    def summarize(
        self, *, query_cell_ids: list[str], probabilities: np.ndarray, labels: list[str]
    ) -> schemas.QueryAnnotationOntologyAwareType:
        """
        Uplift per-cell class probabilities over the cell type ontology.

        :param query_cell_ids: IDs of the query cells, one per row of ``probabilities``.
        :param probabilities: Per-cell class probabilities of shape ``(n_cells, n_labels)``. Each row is
            expected to already be a probability distribution (e.g. a softmax over raw logits).
        :param labels: Class labels corresponding to the probability columns, one per column.

        :raises exceptions.ClassificationLabelSpaceError: If none of ``labels`` maps onto a term in the
            configured cell ontology resource.

        :return: A list of `schemas.QueryCellNeighborhoodOntologyAware`, one per query cell.
        """
        normalized_labels = [normalize_cl_term_id(label) for label in labels]

        if not any(label in self.cell_ontology_resource.ancestors_dictionary for label in normalized_labels):
            raise exceptions.ClassificationLabelSpaceError(
                f"None of the {len(normalized_labels)} model class labels are present in the configured "
                f"cell ontology resource. First labels: {normalized_labels[:5]}"
            )

        result = []
        for i, query_cell_id in enumerate(query_cell_ids):
            scores, unrecognized = accumulate_ontology_scores(
                resource=self.cell_ontology_resource, term_ids=normalized_labels, weights=probabilities[i]
            )
            matches = build_ontology_matches(
                resource=self.cell_ontology_resource, scores=scores, prune_threshold=self.prune_threshold
            )
            result.append(
                schemas.QueryCellNeighborhoodOntologyAware(
                    query_cell_id=query_cell_id,
                    matches=matches,
                    total_weight=float(probabilities[i].sum()),
                    total_neighbors=len(labels),
                    total_neighbors_unrecognized=unrecognized,
                )
            )

        return result
