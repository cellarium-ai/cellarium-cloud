import asyncio
import pickle
import typing as t

import numpy as np
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor
from smart_open import open

from casp.services import settings
from casp.services.api import schemas
from casp.services.api.data_manager import CellOperationsDataManager

DEFAULT_SUMMARIZE_NEIGHBOR_CONTEXT_METHOD = "ontology_aware"


class CellOntologyResource:
    """
    Class to store the cell ontology resource data. This includes:
     - Ancestors dictionary, which has to be precomputed and stored in the bucket
     - Cell ontology term ID to cell type mapping dictionary
    """
    def __init__(self, cell_ontology_resource_dict: t.Dict[str, t.Any]):
        if "ancestors_dictionary" not in cell_ontology_resource_dict:
            raise ValueError("Ancestors dictionary is not found in the cell ontology resource dictionary")
        if "cell_ontology_term_id_to_cell_type" not in cell_ontology_resource_dict:
            raise ValueError(
                "Cell ontology term ID to cell type mapping is not found in the cell ontology resource dictionary"
            )

        self.ancestors_dictionary = cell_ontology_resource_dict["ancestors_dictionary"]
        self.ontology_term_id_to_name_dict = cell_ontology_resource_dict["cell_ontology_term_id_to_cell_type"]


class ConsensusEngine:

    def __init__(self, cell_operations_dm: t.Optional[CellOperationsDataManager] = None):
        self.cell_operations_dm = cell_operations_dm or CellOperationsDataManager()

        with open(settings.GCS_CELL_ONTOLOGY_RESOURCE_FILE, "rb") as f:
            cell_ontology_dict = pickle.load(f)
            self.cell_ontology_resource = CellOntologyResource(cell_ontology_resource_dict=cell_ontology_dict)

    def get_ontology_aware_neighborhood_context(
        self,
        neighbors: t.List[MatchNeighbor],
        neighbors_metadata_dict: t.Dict[str, schemas.CellariumCellMetadata],
        normalize: bool = True,
        prune_threshold: t.Optional[float] = 0.2,
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Get the neighborhood context for a query cell using the ontology-aware method. The ontology-aware method assigns
        weights to each neighbor based on their distance and the cell type ontology.

        :param neighbors: Neighbors of the query cell from the matching engine
        :param neighbors_metadata_dict: Dictionary of neighbor metadata
        :param normalize: Whether to normalize the scores by the total weight
        :param prune_threshold: Threshold to prune the response by. If None, no pruning is done.

        :return: Dictionary of scores for each cell type ontology term
        """
        neighbor_distances = np.asarray([neighbor.distance for neighbor in neighbors])

        neighbor_metadata = [neighbors_metadata_dict[neighbor.id] for neighbor in neighbors]

        # Get weights for each neighbor
        gamma = 1.0 / np.median(neighbor_distances)
        weights = np.exp(gamma * neighbor_distances)

        total_weight = 0
        total_usable_neighbors = 0
        scores_dict = {k: 0 for k in self.cell_ontology_resource.ancestors_dictionary.keys()}

        for neighbor_metadata, weight in zip(neighbor_metadata, weights):
            # Cell Ontology IDs have the formats: CL:0000000 (in our database) and CL_0000000 (in the ontology graph)
            neighbor_cell_type_ontology_id = neighbor_metadata.cell_type_ontology_term_id.replace(":", "_")

            total_weight += weight
            total_usable_neighbors += 1

            # Add weight to the neighbor's ontology ID
            scores_dict[neighbor_cell_type_ontology_id] += weight

            # Propagate the weight to all ancestors (consistent subgraph of the cell type ontology)
            ancestor_ontology_ids = self.cell_ontology_resource.ancestors_dictionary[neighbor_cell_type_ontology_id]

            for ancestor in ancestor_ontology_ids:
                scores_dict[ancestor] += weight

        if normalize:
            scores_dict = {k: v / total_weight for k, v in scores_dict.items()}

        if prune_threshold is not None:
            scores_dict = {k: v for k, v in scores_dict.items() if v >= prune_threshold}

        result = [
            {
                "cell_type_ontology_term_id": k,
                "cell_type": self.cell_ontology_resource.ontology_term_id_to_name_dict[k],
                "score": v,
            }
            for k, v in scores_dict.items()
        ]

        return result

    def summarize_query_neighbor_context(
        self,
        query_ids: t.List[str],
        knn_query: t.List[t.List[MatchNeighbor]],
        method: str = DEFAULT_SUMMARIZE_NEIGHBOR_CONTEXT_METHOD,
        normalize: bool = True,
    ) -> t.List:
        unique_neighbor_ids = set(int(neighbor.id) for query_neighbors in knn_query for neighbor in query_neighbors)
        neighbors_metadata = self.cell_operations_dm.get_cell_metadata_by_ids(list(unique_neighbor_ids))
        neighbors_metadata_dict = {str(neighbor.cas_cell_index): neighbor for neighbor in neighbors_metadata}

        result = []
        for query_id, query_neighbors in zip(query_ids, knn_query):
            match method:
                case "ontology_aware":
                    neighbor_context = self.get_ontology_aware_neighborhood_context(
                        neighbors=query_neighbors,
                        neighbors_metadata_dict=neighbors_metadata_dict,
                        normalize=normalize,
                    )
                case _:
                    raise ValueError(f"Unknown method: {method}")

            result.append({"query_cell_id": query_id, "matches": neighbor_context})

        return result

    async def summarize_query_neighbor_context_async(
        self,
        query_ids: t.List[str],
        knn_query: t.List[t.List[MatchNeighbor]],
        method: str = DEFAULT_SUMMARIZE_NEIGHBOR_CONTEXT_METHOD,
        normalize: bool = True,
    ) -> t.List:
        result = await asyncio.to_thread(
            self.summarize_query_neighbor_context,
            query_ids=query_ids,
            knn_query=knn_query,
            method=method,
            normalize=normalize,
        )
        return result
