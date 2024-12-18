import json
import typing as t

import numpy as np
from smart_open import open

from casp.services import settings
from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.data_manager import CellOperationsDataManager
from casp.services.api.services.consensus_engine.strategies.base import ConsensusStrategyProtocol


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

    def __init__(self, cell_ontology_resource_dict: t.Optional[t.Dict[str, t.Any]] = None):
        if cell_ontology_resource_dict is None:
            with open(settings.GCS_CELL_ONTOLOGY_RESOURCE_FILE, "rb") as f:
                cell_ontology_resource_dict = json.loads(f.read())

        if "ancestors_dictionary" not in cell_ontology_resource_dict:
            raise ValueError("`ancestors_dictionary` is not found in the cell ontology resource dictionary")
        if "cell_ontology_term_id_to_cell_type" not in cell_ontology_resource_dict:
            raise ValueError(
                "`cell_ontology_term_id_to_cell_type` mapping is not found in the cell ontology resource dictionary"
            )

        self.ancestors_dictionary = cell_ontology_resource_dict["ancestors_dictionary"]
        self.ontology_term_id_to_name_dict = cell_ontology_resource_dict["cell_ontology_term_id_to_cell_type"]


class CellTypeOntologyAwareConsensusStrategy(ConsensusStrategyProtocol):
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
    """

    REQUIRED_CELL_INFO_FEATURE_NAMES = ["cas_cell_index", "cell_type", "cell_type_ontology_term_id"]

    def __init__(
        self,
        prune_threshold: float,
        weighting_prefactor: float,
        cell_ontology_resource: t.Optional[CellOntologyResource] = None,
        cell_operations_dm: t.Optional[CellOperationsDataManager] = None,
    ):
        self.cell_ontology_resource = cell_ontology_resource or CellOntologyResource()
        self.cell_operations_dm = cell_operations_dm or CellOperationsDataManager()
        self.prune_threshold = prune_threshold
        self.weighting_prefactor = weighting_prefactor

    def _calculate_cell_type_ontology_aware_scores(
        self,
        query_cell_id: str,
        neighbors: t.List[MatchResult.Neighbor],
        neighbors_metadata_dict: t.Dict[str, schemas.CellariumCellMetadata],
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

        total_weight = 0
        total_neighbors = 0
        total_neighbors_unrecognized = 0
        scores_dict = {k: 0 for k in self.cell_ontology_resource.ancestors_dictionary.keys()}

        for neighbor_metadata, weight in zip(neighbor_metadata, weights):
            # Cell Ontology IDs have the formats: CL:0000000 (in our database) and CL_0000000 (in the ontology graph)
            neighbor_cell_type_ontology_id = neighbor_metadata.cell_type_ontology_term_id.replace(":", "_")

            total_weight += weight
            total_neighbors += 1

            try:
                # Add weight to the neighbor's ontology ID
                scores_dict[neighbor_cell_type_ontology_id] += weight
            except KeyError:
                total_neighbors_unrecognized += 1
                continue
            else:
                ancestor_ontology_ids = self.cell_ontology_resource.ancestors_dictionary[neighbor_cell_type_ontology_id]

                # Propagate the weight to all ancestors (consistent subgraph of the cell type ontology)
                for ancestor in ancestor_ontology_ids:
                    scores_dict[ancestor] += weight

        # Normalize the weights
        scores_dict = {k: v / total_weight for k, v in scores_dict.items()}

        if self.prune_threshold > 0.0:
            scores_dict = {k: v for k, v in scores_dict.items() if v >= self.prune_threshold}

        neighborhood_summary = [
            schemas.NeighborhoodSummaryOntologyAware(
                score=score,
                cell_type_ontology_term_id=cell_type_ontology_term_id,
                cell_type=self.cell_ontology_resource.ontology_term_id_to_name_dict[cell_type_ontology_term_id],
            )
            for cell_type_ontology_term_id, score in scores_dict.items()
        ]
        return schemas.QueryCellNeighborhoodOntologyAware(
            query_cell_id=query_cell_id,
            matches=neighborhood_summary,
            total_weight=total_weight,
            total_neighbors=total_neighbors,
            total_neighbors_unrecognized=total_neighbors_unrecognized,
        )

    def summarize(
        self, query_cell_ids: t.List[str], knn_query: MatchResult
    ) -> schemas.QueryAnnotationOntologyAwareType:
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
            cell_ids=list(unique_neighbor_ids),
            metadata_feature_names=self.REQUIRED_CELL_INFO_FEATURE_NAMES,
        )
        neighbors_metadata_dict = {str(neighbor.cas_cell_index): neighbor for neighbor in neighbors_metadata}

        result = []
        for query_cell_id, query_neighbors in zip(query_cell_ids, knn_query.matches):
            query_cell_neighborhood = self._calculate_cell_type_ontology_aware_scores(
                query_cell_id=query_cell_id,
                neighbors=query_neighbors.neighbors,
                neighbors_metadata_dict=neighbors_metadata_dict,
            )
            result.append(query_cell_neighborhood)

        return result
