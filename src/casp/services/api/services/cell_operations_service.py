"""
Cellarium Service Controller. It provides methods to communicate with services in Cellarium Cloud
infrastructure over different protocols in async manner.
"""

import typing as t

import numpy as np
from google.cloud.aiplatform.matching_engine.matching_engine_index_endpoint import MatchNeighbor

from casp.services.api import clients, schemas
from casp.services.api.data_manager import CellariumGeneralDataManager, CellOperationsDataManager
from casp.services.api.data_manager import exceptions as dm_exc
from casp.services.api.services import exceptions
from casp.services.db import models
from casp.services.utils import numpy_utils

AVAILABLE_FIELDS_DICT = set(schemas.CellariumCellMetadata.__fields__.keys())


class CellOperationsService:
    """
    Cell Analysis Service. It provides methods to communicate with services in Cellarium Cloud infrastructure.
    It leverages async communication with services and uses data access objects to communicate with
    datastore.
    """

    def __init__(self):
        self.cell_operations_dm = CellOperationsDataManager()
        self.cellarium_general_dm = CellariumGeneralDataManager()

    def authorize_model_for_user(self, user: models.User, model_name: str) -> None:
        """
        Authorize user to use a specific model. If user is not authorized, raise an exception.

        :param user: User object to check permissions for.
        :param model_name: Model name to check permissions for.
        """
        if user.is_admin:
            return
        try:
            model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)
        except dm_exc.NotFound as e:
            raise exceptions.InvalidInputError(str(e))
        if model.admin_use_only and not user.is_admin:
            raise exceptions.AccessDeniedError(
                f"{model_name} model is not available. Please reach out to the Cellarium team for more information."
            )

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

    def __get_match_index_endpoint_client_for_model(
        self, model_name: str
    ) -> t.Tuple[models.CASMatchingEngineIndex, clients.CustomMatchingEngineIndexEndpointClient]:
        index = self.cellarium_general_dm.get_index_for_model(model_name=model_name)

        return index, clients.CustomMatchingEngineIndexEndpointClient(index_endpoint_name=index.endpoint_id)

    @staticmethod
    def __validate_knn_response(
        embeddings: np.array, knn_response: t.Union[t.List[t.List[MatchNeighbor]], t.List[t.List[t.Dict[str, t.Any]]]]
    ) -> None:
        if len(knn_response) != len(embeddings):
            raise exceptions.VectorSearchResponseError(
                f"Number of query ids ({len(embeddings)}) and knn matches ({len(knn_response)}) does not match. "
                f"This could probably be caused by Vector Search overload."
            )
        if len(knn_response[0]) == 0:
            raise exceptions.VectorSearchResponseError("Vector Search returned match with 0 neighbors.")

    def get_knn_matches(self, embeddings: np.array, model_name: str) -> t.List[t.List[MatchNeighbor]]:
        """
        Run KNN matching synchronously using Matching Engine client over gRPC.

        :param embeddings: Embeddings from the model.
        :param model_name: Model name to use for matching.

        :return: List of lists of MatchNeighbor objects.
        """
        index, index_endpoint_client = self.__get_match_index_endpoint_client_for_model(model_name=model_name)

        matches = index_endpoint_client.match(
            deployed_index_id=index.deployed_index_id,
            queries=embeddings,
            num_neighbors=index.num_neighbors,
        )
        self.__validate_knn_response(embeddings=embeddings, knn_response=matches)
        return matches

    def get_knn_matches_as_dict(self, embeddings: np.array, model_name: str) -> t.List[t.List[t.Dict[str, t.Any]]]:
        index, index_endpoint_client = self.__get_match_index_endpoint_client_for_model(model_name=model_name)

        matches = index_endpoint_client.match_as_dict(
            deployed_index_id=index.deployed_index_id,
            queries=embeddings,
            num_neighbors=index.num_neighbors,
        )
        self.__validate_knn_response(embeddings=embeddings, knn_response=matches)
        return matches

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

        temp_table_fqn = self.cell_operations_dm.insert_matches_to_temp_table(
            query_ids=query_ids, knn_response=knn_response
        )

        if include_dev_metadata:
            return self.cell_operations_dm.get_neighborhood_distance_summary_dev_details(
                cas_model=cas_model, match_temp_table_fqn=temp_table_fqn
            )

        return self.cell_operations_dm.get_neighborhood_distance_summary(
            cas_model=cas_model, match_temp_table_fqn=temp_table_fqn
        )

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
        self.authorize_model_for_user(user=user, model_name=model_name)
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

    async def search_adata_file(
        self, user: models.User, file: t.BinaryIO, model_name: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Search for similar cells in a single anndata file with Cellarium CAS. Input file should be validated and
        sanitized according to the model schema.

        :param user: User object to check permissions for.
        :param file: Byte object of :class:`anndata.AnnData` file to search.
        :param model_name: Model name to use for search. See `/list-models` endpoint for available models.

        :return: JSON response with search results.
        """
        self.authorize_model_for_user(user=user, model_name=model_name)
        query_ids, embeddings = await self.get_embeddings(file_to_embed=file, model_name=model_name)
        knn_matches_dict = self.get_knn_matches_as_dict(embeddings=embeddings, model_name=model_name)

        return [
            {
                "query_cell_id": query_ids[i],
                "neighbors": knn_matches_dict[i],
            }
            for i in range(0, len(query_ids))
        ]

    def get_cells_by_ids_for_user(
        self, user: models.User, cell_ids: t.List[int], metadata_feature_names: t.List[str], model_name: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Get cells by their ids from BigQuery `cas_cell_info` table.

        :param user: User object to check permissions for
        :param cell_ids: Cas cell indexes from BigQuery
        :param metadata_feature_names: Metadata features to return from BigQuery `cas_cell_info` table
        :param model_name: Name of the model to query. Used to get the dataset name where to get the cells from

        :return: List of dictionaries representing the query results.
        """
        self.authorize_model_for_user(user=user, model_name=model_name)
        for feature_name in metadata_feature_names:
            if feature_name not in AVAILABLE_FIELDS_DICT:
                raise exceptions.CellMetadataColumnDoesntExist(f"Feature {feature_name} is not available for querying")

        if "cas_cell_index" not in metadata_feature_names:
            metadata_feature_names.append("cas_cell_index")
        try:
            return self.cell_operations_dm.get_cell_metadata_by_ids(
                cell_ids=cell_ids,
                metadata_feature_names=metadata_feature_names,
                model_name=model_name,
            )
        except dm_exc.NotFound as e:
            raise exceptions.InvalidInputError(str(e))
