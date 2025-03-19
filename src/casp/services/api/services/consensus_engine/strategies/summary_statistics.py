import typing as t
from collections import OrderedDict, defaultdict

import numpy as np

from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchResult
from casp.services.api.data_manager import CellOperationsDataManager
from casp.services.api.services.consensus_engine.strategies.base import ConsensusStrategyInterface


class CellTypeSummaryStatisticsConsensusStrategy(ConsensusStrategyInterface):
    """
    Handle cell type count consensus strategy, summarizing query neighbor context by cell type distribution.
    """

    REQUIRED_CELL_INFO_FEATURE_NAMES = ["cas_cell_index", "cell_type"]

    def __init__(self, cell_operations_dm: CellOperationsDataManager):
        self.cell_operations_dm = cell_operations_dm

    @staticmethod
    def calculate_cell_type_summary_stats_for_query_cell(
        query_cell_id: str,
        neighbors: t.List[MatchResult.Neighbor],
        neighbors_metadata_dict: t.Dict[str, schemas.CellariumCellMetadata],
    ) -> schemas.QueryCellNeighborhoodCellTypeSummaryStatistics:
        """
        Calculate summary statistics for each cell_type based on nearest neighbor search

        :param query_cell_id: Querying cell ID
        :param neighbors: A list of neighbors from nearest neighbor search
        :param neighbors_metadata_dict:  A dictionary with O1 run-time access to metadata by cas_cell_index

        :return: Summary Statistics for each cell type based on nearest neighbor search
        """
        grouped_distances = defaultdict(list)

        # Group distances by cell_type in-place
        for neighbor in neighbors:
            cell_type = neighbors_metadata_dict[neighbor.cas_cell_index].cell_type
            grouped_distances[cell_type].append(neighbor.distance)

        # Order grouped_distances by the number of matches (length of distances)
        ordered_grouped_distances = OrderedDict(
            sorted(grouped_distances.items(), key=lambda x: len(x[1]), reverse=True)
        )
        query_cell_matches = []

        for cell_type, distances in ordered_grouped_distances.items():
            distances_np = np.array(distances)
            query_cell_matches.append(
                schemas.NeighborhoodCellTypeSummaryStatistics(
                    cell_type=cell_type,
                    cell_count=len(distances_np),
                    min_distance=distances_np.min(),
                    p25_distance=np.percentile(a=distances_np, q=25),
                    median_distance=np.median(distances_np),
                    p75_distance=np.percentile(a=distances_np, q=75),
                    max_distance=distances_np.max(),
                )
            )
        return schemas.QueryCellNeighborhoodCellTypeSummaryStatistics(
            query_cell_id=query_cell_id, matches=query_cell_matches
        )

    def summarize(
        self, query_cell_ids: t.List[str], knn_response: MatchResult
    ) -> t.List[schemas.QueryCellNeighborhoodCellTypeSummaryStatistics]:
        """
        Summarize query neighbor context by cell type distribution, querying a database for the distribution.

        :param query_cell_ids: List of query cell IDs.
        :param knn_response: The result of the kNN query.

        :return: An instance of `schemas.QueryAnnotationCellTypeCountType`, representing the summarized cell type
            distribution of query neighbors.
        """
        unique_neighbor_ids = knn_response.get_unique_ids()
        neighbors_metadata = self.cell_operations_dm.get_cell_metadata_by_ids(
            cell_ids=unique_neighbor_ids,
            metadata_feature_names=self.REQUIRED_CELL_INFO_FEATURE_NAMES,
        )
        neighbors_metadata_dict = {str(neighbor.cas_cell_index): neighbor for neighbor in neighbors_metadata}

        results = []
        for query_cell_id, query_neighbors in zip(query_cell_ids, knn_response.matches):
            query_cell_neighborhood = self.calculate_cell_type_summary_stats_for_query_cell(
                query_cell_id=query_cell_id,
                neighbors=query_neighbors.neighbors,
                neighbors_metadata_dict=neighbors_metadata_dict,
            )
            results.append(query_cell_neighborhood)

        return results
