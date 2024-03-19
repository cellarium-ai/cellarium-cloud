import typing as t
import numpy as np
import scipy.sparse as sp
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor
from casp.services.api.data_manager import CellOperationsDataManager
from casp.services.api import schemas
from casp.services import settings
from smart_open import open

DEFAULT_SUMMARIZE_NEIGHBOR_CONTEXT_METHOD = "ontology_aware"


class ConsensusEngine:

    def __init__(self):
        self.cell_operations_dm = CellOperationsDataManager()

        with open(settings.GCS_CELL_ONTOLOGY_SPARSE_MATRIX_NPZ, "rb") as f:
            # self.cell_ontology_csr_matrix = sp.load_npz(f)
        self.cell_ontology_csr_matrix = sp.load_npz("cell_ontology_sparse_matrix.npz")

    # def summarize_neighbor_context_ontology_aware_method(
    #     self, query_ids: t.List[str], query_neighbors: t.List[t.List[MatchNeighbor]]
    # ) -> t.List:
    #     unique_neighbor_ids = set()
    #
    #     for neighbors in query_neighbors:
    #         for neighbor in neighbors:
    #             unique_neighbor_ids.add(neighbor.id)

    def ontology_aware_summarize_neighbor_context(
        self, neighbors: t.List[MatchNeighbor], neighbor_metadata_dict: t.Dict[str, schemas.CellariumCellMetadata]
    ) -> t.Dict[str, t.Any]:
        neighbor_distances = np.asarray([neighbor.distance for neighbor in neighbors])
        neighbor_metadata = [neighbor_metadata_dict[neighbor.id] for neighbor in neighbors]
        gamma = 1. / np.median(neighbor_distances)
        weights = np.exp(gamma * neighbor_distances)

        total_weight = 0
        total_usable_neighbors = 0
        # scores_array =

        for weight, neighbor_metadata in zip(weights, neighbor_metadata):
            neighbor_metadata["weight"] = weight

            neighbor_cl_name = neighbor_metadata.cell_type_ontology_term_id.replace(':', '_')

            if neighbor_cl_name not in allowed_cl_names_set:
                continue

            nn_cl_idx = cl_names_to_idx_map[nn_cl_name]
            total_usable_neighbors += 1
            total_weight += weight
            scores_array[nn_cl_idx] += weight
            # propagate the weight to all ancestors (consistent subgraph of the CL onotlogy)
            scores_array[cl_ancestors_csr_matrix[nn_cl_idx].indices] += weight



    def summarize_query_neighbor_context(
        self,
        query_ids: t.List[str],
        knn_query: t.List[t.List[MatchNeighbor]],
        method: str = DEFAULT_SUMMARIZE_NEIGHBOR_CONTEXT_METHOD
    ) -> t.List:
        unique_neighbor_ids = set(int(neighbor.id) for query_neighbors in knn_query for neighbor in query_neighbors)
        neighbor_metadata = self.cell_operations_dm.get_cell_metadata_by_ids(list(unique_neighbor_ids))
        neighbor_metadata_dict = {neighbor.cas_cell_index: neighbor for neighbor in neighbor_metadata}

        result = []
        for query_id, query_neighbors in zip(query_ids, knn_query):
            match method:
                case "ontology_aware":
                    neighbor_context = self.ontology_aware_summarize_neighbor_context(query_neighbors)
                case _:
                    raise ValueError(f"Unknown method: {method}")

            result.append({"query_cell_id": query_id, "matches": neighbor_context})

        return result
