"""
Cellarium Service Controller. It provides methods to communicate with services in Cellarium Cloud
infrastructure over different protocols in async manner.
"""
import typing as t

import numpy as np
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor

from casp.services.api import clients
from casp.services.api.data_manager import CellAnalysisDataManager, CellariumGeneralDataManager
from casp.services.db import models
from casp.services.utils import numpy_utils


class CellAnalysisService:
    """
    Cell Analysis Service. It provides methods to communicate with services in Cellarium Cloud infrastructure.
    It leverages async communication with services and uses data access objects to communicate with
    datastore.
    """

    def __init__(self):
        self.cell_analysis_dm = CellAnalysisDataManager()
        self.cellarium_general_dm = CellariumGeneralDataManager()

    @staticmethod
    async def get_embeddings(file_to_embed: t.BinaryIO, model_name: str) -> t.Tuple[t.List[str], np.array]:
        """
        Get embeddings from model inference service. Unwrap response and return query ids and embeddings.
        Since model embedding service returns embeddings in base64 format, we need to convert it to numpy array.

        :param file_to_embed: File object of :class:`anndata.AnnData` object to embed.
        :param model_name: Model name to use for embedding.

        :return: Query ids (original cell ids from the input file) and embeddings.
        """
        embeddings_response_json = await clients.ModelInferenceClient.call_model_embed(
            file_to_embed=file_to_embed.read(), model_name=model_name
        )
        query_ids = embeddings_response_json["obs_ids"]
        embeddings = numpy_utils.base64_to_numpy(embeddings_response_json["embeddings_b64"])
        return query_ids, embeddings

    def get_knn_matches(self, embeddings: np.array, model_name: str) -> t.List[t.List[MatchNeighbor]]:
        """
        Run KNN matching synchronously using Matching Engine client over gRPC.

        :param embeddings: Embeddings from the model.
        :param model_name: Model name to use for matching.

        :return: List of lists of MatchNeighbor objects.
        """
        cas_model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)
        cas_matching_engine_index = cas_model.cas_matching_engine

        index_endpoint = clients.CustomMatchingEngineIndexEndpointClient(
            index_endpoint_name=cas_matching_engine_index.endpoint_id
        )

        return index_endpoint.match(
            deployed_index_id=cas_matching_engine_index.deployed_index_id,
            queries=embeddings,
            num_neighbors=cas_matching_engine_index.num_neighbors,
        )

    def get_cell_type_distribution(
        self,
        query_ids: t.List[str],
        knn_response: t.List[t.List[MatchNeighbor]],
        model_name: str,
        include_dev_metadata: bool,
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Get cell type distribution for the given query ids and knn response. Put results into temporary table and
        aggregate summary statistics.

        :param query_ids: List of query ids (original cell ids from the input file).
        :param knn_response: List of lists of MatchNeighbor objects returned by space vector search service.
        :param model_name: Model name to use for matching.
        :param include_dev_metadata: Boolean flag to include dev metadata in the response.

        :return: List of cell type distribution objects.
        """
        cas_model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)

        temp_table_fqn = self.cell_analysis_dm.insert_matches_to_temp_table(
            query_ids=query_ids, knn_response=knn_response
        )

        if include_dev_metadata:
            return self.cell_analysis_dm.get_match_query_metadata_dev_details(
                cas_model=cas_model, match_temp_table_fqn=temp_table_fqn
            )

        return self.cell_analysis_dm.get_match_query_metadata(cas_model=cas_model, match_temp_table_fqn=temp_table_fqn)

    async def annotate_adata_file(
        self, user: models.User, file: t.BinaryIO, model_name: str, include_dev_metadata: bool
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
        :param include_dev_metadata: Boolean flag indicating whether to include dev metadata in the response.

        :return: JSON response with annotations.
        """
        query_ids, embeddings = await self.get_embeddings(file_to_embed=file, model_name=model_name)
        knn_response = self.get_knn_matches(embeddings=embeddings, model_name=model_name)
        annotation_response = self.get_cell_type_distribution(
            query_ids=query_ids,
            knn_response=knn_response,
            model_name=model_name,
            include_dev_metadata=include_dev_metadata,
        )
        self.cellarium_general_dm.increment_user_cells_processed(user=user, number_of_cells=len(query_ids))
        return annotation_response
