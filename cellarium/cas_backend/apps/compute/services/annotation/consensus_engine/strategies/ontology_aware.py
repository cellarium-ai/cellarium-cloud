import numpy as np

from cellarium.cas_backend.apps.compute import schemas
from cellarium.cas_backend.apps.compute.services.annotation.consensus_engine.strategies.base import (
    ConsensusStrategyInterface,
)
from cellarium.cas_backend.apps.compute.services.annotation.ontology import (
    CellOntologyResource,
    accumulate_ontology_scores,
    build_ontology_matches,
)
from cellarium.cas_backend.apps.compute.vector_search import MatchResult
from cellarium.cas_backend.core.data_managers import CellOperationsDataManager


class CellTypeOntologyAwareConsensusStrategy(ConsensusStrategyInterface):
    """
    Handle ontology-aware consensus strategy, summarizing query neighbor context by cell type ontology. Weights are
    assigned to each neighbor cell type based on their distance and cell type ontology, and the weights are propagated
    to their ancestors in the cell type ontology graph.

    Algorithm:

    1. Get metadata for each unique neighbor cell.
    2. Iterate over each query cell.
       2.1. Get distances for each neighbor cell.
       2.2. Calculate weights for each neighbor cell.
       2.3. Update weights for each neighbor cell and their ancestors in the cell type ontology.
       2.4. Normalize the weights

    :param prune_threshold: Threshold for pruning weights below a certain value. If 0, no pruning is performed.
    :param cell_ontology_resource: Cell ontology resource object.
    :param cell_operations_dm: Cell operations data manager object.
    :param weighting_prefactor: Distance exponential weighting prefactor.
    :param cell_metadata_uri: GCS URI pointing to the TileDB SOMA DataFrame for this model's cell metadata.
    """

    REQUIRED_CELL_INFO_FEATURE_NAMES = ["cas_cell_index", "cell_type", "cell_type_ontology_term_id"]

    def __init__(
        self,
        prune_threshold: float,
        weighting_prefactor: float,
        cell_ontology_resource: CellOntologyResource,
        cell_metadata_uri: str,
        cell_operations_dm: CellOperationsDataManager | None = None,
    ):
        self.cell_ontology_resource = cell_ontology_resource
        self.cell_operations_dm = cell_operations_dm or CellOperationsDataManager()
        self.cell_metadata_uri = cell_metadata_uri
        self.prune_threshold = prune_threshold
        self.weighting_prefactor = weighting_prefactor

    def _calculate_cell_type_ontology_aware_scores(
        self,
        query_cell_id: str,
        neighbors: list[MatchResult.Neighbor],
        neighbors_metadata_dict: dict[str, schemas.CellariumCellMetadata],
    ) -> schemas.QueryCellNeighborhoodOntologyAware:
        """
        Utilize the ontology-aware method to assign weights to neighbor cells based on their distance and cell type
        ontology, to inform context summarization.

        :param query_cell_id: ID of the query cell.
        :param neighbors: Neighbors of the query cell as determined by the matching engine.
        :param neighbors_metadata_dict: Metadata for each neighbor cell.
        :return: A list of `schemas.AnnotationInfoOntologyAware` instances, each representing the weighted cell type
             ontology term for a neighbor cell.
        """

        neighbor_distances = np.asarray([neighbor.distance for neighbor in neighbors])

        neighbor_metadata = [neighbors_metadata_dict[neighbor.cas_cell_index] for neighbor in neighbors]

        # Get weights for each neighbor
        gamma = -self.weighting_prefactor / np.median(neighbor_distances)
        weights = np.exp(gamma * neighbor_distances)

        term_ids = [metadata.cell_type_ontology_term_id for metadata in neighbor_metadata]
        scores, total_neighbors_unrecognized = accumulate_ontology_scores(
            resource=self.cell_ontology_resource, term_ids=term_ids, weights=weights
        )

        total_weight = float(weights.sum())
        total_neighbors = len(neighbors)

        # Normalize the weights
        scores = {k: v / total_weight for k, v in scores.items()}

        matches = build_ontology_matches(
            resource=self.cell_ontology_resource, scores=scores, prune_threshold=self.prune_threshold
        )
        return schemas.QueryCellNeighborhoodOntologyAware(
            query_cell_id=query_cell_id,
            matches=matches,
            total_weight=total_weight,
            total_neighbors=total_neighbors,
            total_neighbors_unrecognized=total_neighbors_unrecognized,
        )

    def summarize(self, query_cell_ids: list[str], knn_query: MatchResult) -> schemas.QueryAnnotationOntologyAwareType:
        """
        Summarize the query neighbor context using the ontology-aware method, assigning weights to each neighbor cell
        type and propagating the weights to their ancestors in the cell type ontology.

        :param query_cell_ids: IDs of the query cells.
        :param knn_query: The result of the kNN query.

        :return: A list of `schemas.QueryCellAnnotationOntologyAware`, representing the summarized context for each
            query cell.
        """
        unique_neighbor_ids = knn_query.get_unique_ids()
        neighbors_metadata = self.cell_operations_dm.get_cell_metadata_by_ids(
            cell_metadata_uri=self.cell_metadata_uri,
            cell_ids=list(unique_neighbor_ids),
            metadata_feature_names=self.REQUIRED_CELL_INFO_FEATURE_NAMES,
        )
        neighbors_metadata_dict = {neighbor.cas_cell_index: neighbor for neighbor in neighbors_metadata}

        result = []
        for query_cell_id, query_neighbors in zip(query_cell_ids, knn_query.matches, strict=False):
            query_cell_neighborhood = self._calculate_cell_type_ontology_aware_scores(
                query_cell_id=query_cell_id,
                neighbors=query_neighbors.neighbors,
                neighbors_metadata_dict=neighbors_metadata_dict,
            )
            result.append(query_cell_neighborhood)

        return result
