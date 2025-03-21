"""
Cellarium Service Controller. It provides methods to communicate with services in Cellarium Cloud
infrastructure over different protocols in async manner.
"""

import logging
import typing as t

import anndata
import numpy as np
from tenacity import AsyncRetrying, stop_after_attempt, wait_exponential, wait_random

from casp.services import settings
from casp.services.api import schemas
from casp.services.api.clients.matching_client import MatchingClient, MatchResult
from casp.services.api.data_manager import CellariumGeneralDataManager, CellOperationsDataManager
from casp.services.api.data_manager import exceptions as dm_exc
from casp.services.api.services import Authorizer, consensus_engine, exceptions
from casp.services.api.services.cell_quota_service import CellQuotaService
from casp.services.db import models
from casp.services.model_inference import services

AVAILABLE_FIELDS_DICT = set(schemas.CellariumCellMetadata.__fields__.keys())

logger = logging.getLogger(__name__)


class CellOperationsService:
    """
    Cell Analysis Service. It provides methods to communicate with services in Cellarium Cloud infrastructure.
    It leverages async communication with services and uses data access objects to communicate with
    datastore.
    """

    def __init__(
        self,
        cell_operations_dm: t.Optional[CellOperationsDataManager] = None,
        cellarium_general_dm: t.Optional[CellariumGeneralDataManager] = None,
        cell_quota_service: t.Optional[CellQuotaService] = None,
        model_service: t.Optional[services.ModelInferenceService] = None,
        authorizer: t.Optional[Authorizer] = None,
    ):
        self.cell_operations_dm = cell_operations_dm or CellOperationsDataManager()
        self.cellarium_general_dm = cellarium_general_dm or CellariumGeneralDataManager()
        self.cell_quota_service = cell_quota_service or CellQuotaService()
        self.model_service = model_service or services.ModelInferenceService()
        self.authorizer = authorizer or Authorizer(cellarium_general_dm=self.cellarium_general_dm)

    @staticmethod
    def get_cache_info(user: models.User) -> schemas.CacheInfo:
        """
        Returns the cache info for the file and module cache in ModelInferenceService.

        :return: A dict containing two entries: file_cache_info and module_cache_info.  Each
            contains a tuple with 4 values: hits, misses, maxsize, and currsize.
        """
        if not user.is_admin:
            raise exceptions.AccessDeniedError("Access denied.  This is an admin-only endpoint.")

        cache_info = services.ModelInferenceService.get_cache_info()

        return schemas.CacheInfo(
            file_cache_info=cache_info["file_cache_info"],
            module_cache_info=cache_info["module_cache_info"],
        )

    def __validate_anndata_file(self, adata: anndata.AnnData):
        """
        Validates the type of the anndata matrix and the presence of the 'total_mrna_umis' column,
        which are required for the data to be processed properly by our cellarium-ml models.

        :param adata: Anndata object to validate

        :raises exceptions.InvalidAnndataFileException: If the anndata matrix is not float32
        """
        if adata.X.dtype != np.float32:
            raise exceptions.InvalidAnndataFileException(
                f"CAS only support anndata matrices with dtype 'float32'. The provided matrix has type '{adata.X.dtype}'."
            )
        if "total_mrna_umis" not in adata.obs:
            raise exceptions.InvalidAnndataFileException(
                "Provided matrix does not contain 'total_mrna_umis' column. This column is required for CAS to work."
            )
        elif adata.obs["total_mrna_umis"].dtype != np.float32:
            raise exceptions.InvalidAnndataFileException(
                f"Provided matrix column 'total_mrna_umis' should be type 'float32', but has type '{adata.obs['total_mrna_umis'].dtype}'."
            )

    def __read_and_validate_anndata_file(self, file: t.BinaryIO) -> anndata.AnnData:
        """
        Reads an anndata file and validates that the matrix is float32

        :param file: Anndata file to read

        :return: Anndata object
        """
        adata = anndata.read_h5ad(file)
        self.__validate_anndata_file(adata)

        # Reset the file pointer to the beginning of the file because we're gonna need to read it
        # again later.
        file.seek(0)

        return adata

    def __verify_quota_and_log_activity(
        self, user: models.User, adata: anndata.AnnData, model_name: str, method: str
    ) -> int:
        """
        Verify that the number of cells in the anndata object does not exceed the user's remaining
        quota and log the user activity in the database.

        :param user: User object to check quota for.
        :param adata: Anndata object to check the number of cells in.
        :param model_name: Model name to log in the database.
        :param method: Method name to log in the database.

        :return: Number of cells in the anndata object.

        :raises exceptions.QuotaExceededException: If the number of cells in the anndata object exceeds the user's quota.
        """
        cell_count = len(adata)
        user_quota = self.cell_quota_service.get_remaining_quota_for_user(user=user)
        if cell_count > user_quota:
            raise exceptions.QuotaExceededException(
                f"User quota exceeded. You have {user_quota} cells left, "
                f"but you are trying to process {cell_count} cells."
            )
        self.cellarium_general_dm.log_user_activity(
            user_id=user.id,
            model_name=model_name,
            method=method,
            cell_count=cell_count,
            event=models.UserActivityEvent.STARTED,
        )

        return cell_count

    async def get_embeddings(self, adata: anndata.AnnData, model: models.CASModel) -> t.Tuple[t.List[str], np.array]:
        """
        Get embeddings from the specified model.

        :param adata: Object of :class:`anndata.AnnData` to embed.
        :param model: Model to use for embedding.

        :return: Query ids (original cell ids from the input anndata) and embeddings.
        """

        query_ids, embeddings = self.model_service.embed_adata(adata=adata, model=model)
        return query_ids, embeddings

    @staticmethod
    def __validate_knn_response(embeddings: np.array, knn_response: MatchResult) -> None:
        if len(knn_response.matches) != len(embeddings):
            raise exceptions.VectorSearchResponseError(
                f"Number of query ids ({len(embeddings)}) and knn matches ({len(knn_response.matches)}) does not match. "
                f"This could probably be caused by Vector Search overload."
            )

        for match in knn_response.matches:
            if len(match.neighbors) == 0:
                raise exceptions.VectorSearchResponseError("Vector Search returned a match with 0 neighbors.")

    async def __get_knn_matches_from_embeddings_with_retry(
        self, embeddings: np.array, model: models.CASModel, chunk_matches_function: t.Callable
    ) -> MatchResult:
        """
        Get KNN matches for embeddings broken into chunks with retry logic

        :param embeddings: Embeddings to match.
        :param model: Model to use for matching.

        :return: List of lists of MatchNeighbor objects.
        """
        matching_client = MatchingClient.from_index(model.cas_matching_engine)

        # Break embeddings into chunks so we don't overload the matching engine
        embeddings_chunks = CellOperationsService.__split_embeddings_into_chunks(
            embeddings=embeddings, chunk_size=settings.GET_MATCHES_CHUNK_SIZE
        )

        all_matches = MatchResult()
        for i in range(0, len(embeddings_chunks)):
            # Set up retry logic for the matching engine requests
            retryer = AsyncRetrying(
                stop=stop_after_attempt(settings.GET_MATCHES_MAX_RETRIES),
                wait=wait_exponential(
                    multiplier=settings.GET_MATCHES_RETRY_BACKOFF_MULTIPLIER,
                    min=settings.GET_MATCHES_RETRY_BACKOFF_MIN,
                    max=settings.GET_MATCHES_RETRY_BACKOFF_MAX,
                )
                + wait_random(0, 2),
                reraise=True,
            )
            matches = await retryer(
                chunk_matches_function, embeddings_chunk=embeddings_chunks[i], client=matching_client
            )
            all_matches = all_matches.concat(matches)

        return all_matches

    @staticmethod
    async def __get_knn_matches_for_chunk(embeddings_chunk: np.array, client: MatchingClient) -> MatchResult:
        """
        Get KNN matches for a chunk of embeddings (split out from get_knn_matches, so it can be
        retried and called in chunks).

        :param embeddings_chunk: Chunk of embeddings to match.
        :param client: Matching engine client that will handle communicating with the matching engine.

        :return: List of lists of MatchNeighbor objects.
        """
        matches = await client.match(queries=embeddings_chunk)
        CellOperationsService.__validate_knn_response(embeddings=embeddings_chunk, knn_response=matches)
        return matches

    async def get_knn_matches_from_embeddings(self, embeddings: np.array, model: models.CASModel) -> MatchResult:
        """
        Run KNN matching synchronously using Matching Engine client over gRPC.

        :param embeddings: Embeddings from the model.
        :param model: Model to use for matching.

        :return: Matches in a MatchResult object
        """

        return await self.__get_knn_matches_from_embeddings_with_retry(
            embeddings=embeddings,
            model=model,
            chunk_matches_function=self.__get_knn_matches_for_chunk,
        )

    @staticmethod
    def __split_embeddings_into_chunks(embeddings: np.array, chunk_size: int) -> t.List[np.array]:
        """
        Splits the embeddings into chunks based on the configured chunk size.

        :param embeddings: The embeddings to split into chunks.
        :param chunk_size: The number of query embeddings to include in each chunk.

        :return: A list of numpy arrays, each containing a chunk of the embeddings.
        """
        if chunk_size <= 0:
            raise ValueError("chunk_size must be greater than 0.")
        if len(embeddings) == 0:
            return []

        return [embeddings[i : i + chunk_size] for i in range(0, len(embeddings), chunk_size)]

    async def get_knn_matches(
        self, adata: anndata.AnnData, model: models.CASModel
    ) -> t.Tuple[t.List[str], MatchResult]:
        """
        Get KNN matches for anndata object.

        :param adata: Anndata object to get KNN matches for.
        :param model: Model to use for matching.

        :return: Matches in a MatchResult object
        """
        logger.info("Getting embeddings")
        query_ids, embeddings = await self.get_embeddings(adata=adata, model=model)

        if embeddings.size == 0:
            # No further processing needed if there are no embeddings, return empty match result
            return [], MatchResult()

        logger.info("Getting KNN matches")
        knn_response = await self.get_knn_matches_from_embeddings(embeddings=embeddings, model=model)

        return query_ids, knn_response

    async def annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        self,
        user: models.User,
        file: t.BinaryIO,
        model_name: str,
    ) -> schemas.QueryAnnotationCellTypeSummaryStatisticsType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.  This is a wrapper method
        for annotate_cell_type_summary_statistics_strategy that logs to the user activity table and verifies the input
        file doesn't exceed the user's remaining quota

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.

        :return: JSON response with annotations.
        """
        model = self.authorizer.authorize_model_for_user(user=user, model_name=model_name)

        # Read the anndata file and validate it
        adata = self.__read_and_validate_anndata_file(file=file)

        # Make sure the user has enough quota to process the file
        cell_count = self.__verify_quota_and_log_activity(
            user=user, adata=adata, model_name=model_name, method="annotate_cell_type_summary_statistics_strategy"
        )

        # Annotate the file and log successful activity (or log failure if an exception is raised)
        try:
            annotation_response = await self.__annotate_cell_type_summary_statistics_strategy(
                adata=adata,
                model=model,
            )

            self.cellarium_general_dm.log_user_activity(
                user_id=user.id,
                model_name=model_name,
                method="annotate_cell_type_summary_statistics_strategy",
                cell_count=cell_count,
                event=models.UserActivityEvent.SUCCEEDED,
            )

            return annotation_response
        except Exception as e:
            try:
                self.cellarium_general_dm.log_user_activity(
                    user_id=user.id,
                    model_name=model_name,
                    method="annotate_cell_type_summary_statistics_strategy",
                    cell_count=cell_count,
                    event=models.UserActivityEvent.FAILED,
                )
            except Exception as e2:
                logger.error(f"Failed to log user activity: {e2}")
            raise e

    async def __annotate_cell_type_summary_statistics_strategy(
        self,
        adata: anndata.AnnData,
        model: models.CASModel,
    ) -> schemas.QueryAnnotationCellTypeSummaryStatisticsType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.

        :param adata: Instance of :class:`anndata.AnnData` to annotate.
        :param model: Model to use for annotation. See `/list-models` endpoint for available models.

        :return: JSON response with annotations.
        """
        query_ids, knn_response = await self.get_knn_matches(adata=adata, model=model)

        strategy = consensus_engine.CellTypeSummaryStatisticsConsensusStrategy(
            cell_operations_dm=self.cell_operations_dm,
        )

        engine = consensus_engine.ConsensusEngine(strategy=strategy)

        return engine.summarize(query_ids=query_ids, knn_query=knn_response)

    async def annotate_cell_type_ontology_aware_strategy_with_activity_logging(
        self, user: models.User, file: t.BinaryIO, model_name: str, prune_threshold: float, weighting_prefactor: float
    ) -> schemas.QueryAnnotationOntologyAwareType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation. This is a wrapper method
        for annotate_cell_type_ontology_aware_strategy that logs to the user activity table and verifies the input
        file doesn't exceed the user's remaining quota

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
        :param prune_threshold: Prune threshold for the ontology-aware annotation strategy.
        :param weighting_prefactor: Distance exponential weighting prefactor.

        :return: JSON response with annotations.

        """
        model = self.authorizer.authorize_model_for_user(user=user, model_name=model_name)

        # Read the anndata file and validate it
        adata = self.__read_and_validate_anndata_file(file=file)

        # Make sure the user has enough quota to process the file
        cell_count = self.__verify_quota_and_log_activity(
            user=user, adata=adata, model_name=model_name, method="annotate_cell_type_ontology_aware_strategy"
        )

        # Annotate the file and log successful activity (or log failure if an exception is raised)
        try:
            annotation_response = await self.__annotate_cell_type_ontology_aware_strategy(
                adata=adata,
                model=model,
                prune_threshold=prune_threshold,
                weighting_prefactor=weighting_prefactor,
            )

            self.cellarium_general_dm.log_user_activity(
                user_id=user.id,
                model_name=model_name,
                method="annotate_cell_type_ontology_aware_strategy",
                cell_count=cell_count,
                event=models.UserActivityEvent.SUCCEEDED,
            )

            return annotation_response
        except Exception as e:
            try:
                self.cellarium_general_dm.log_user_activity(
                    user_id=user.id,
                    model_name=model_name,
                    method="annotate_cell_type_ontology_aware_strategy",
                    cell_count=cell_count,
                    event=models.UserActivityEvent.FAILED,
                )
            except Exception as e2:
                logger.error(f"Failed to log user activity: {e2}")
            raise e

    async def __annotate_cell_type_ontology_aware_strategy(
        self, adata: anndata.AnnData, model: models.CASModel, prune_threshold: float, weighting_prefactor: float
    ) -> schemas.QueryAnnotationOntologyAwareType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.

        :param user: User object used to increment user cells processed counter.
        :param adata: Object of :class:`anndata.AnnData` to annotate.
        :param model: Model to use for annotation. See `/list-models` endpoint for available models.
        :param prune_threshold: Prune threshold for the ontology-aware annotation strategy.
        :param weighting_prefactor: Distance exponential weighting prefactor.

        :return: JSON response with annotations.
        """
        query_ids, knn_response = await self.get_knn_matches(adata=adata, model=model)

        logger.info("Applying CellTypeOntologyAwareConsensusStrategy to the query results")
        strategy = consensus_engine.CellTypeOntologyAwareConsensusStrategy(
            cell_operations_dm=self.cell_operations_dm,
            prune_threshold=prune_threshold,
            weighting_prefactor=weighting_prefactor,
        )

        logger.info("Summarizing query neighbor context using the specified strategy")
        engine = consensus_engine.ConsensusEngine(strategy=strategy)

        logger.info("Performing final summarization")
        try:
            return engine.summarize(query_ids=query_ids, knn_query=knn_response)
        except dm_exc.CellMetadataDatabaseError as e:
            raise exceptions.InvalidInputError(str(e))

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
        # Read the anndata file and validate it
        adata = self.__read_and_validate_anndata_file(file=file)

        model = self.authorizer.authorize_model_for_user(user=user, model_name=model_name)
        query_ids, embeddings = await self.get_embeddings(adata=adata, model=model)
        if embeddings.size == 0:
            # No further processing needed if there are no embeddings
            return []
        knn_response = await self.get_knn_matches_from_embeddings(embeddings=embeddings, model=model)

        return [
            {
                "query_cell_id": query_ids[i],
                "neighbors": [
                    {
                        "cas_cell_index": neighbor.cas_cell_index,
                        "distance": neighbor.distance,
                    }
                    for neighbor in knn_response.matches[i].neighbors
                ],
            }
            for i in range(0, len(query_ids))
        ]

    def get_cells_by_ids_for_user(
        self, cell_ids: t.List[int], metadata_feature_names: t.List[str]
    ) -> t.List[schemas.CellariumCellMetadata]:
        """
        Get cells by their ids from BigQuery `cas_cell_info` table.

        :param user: User object to check permissions for
        :param cell_ids: Cas cell indexes from BigQuery
        :param metadata_feature_names: Metadata features to return from BigQuery `cas_cell_info` table

        :return: List of dictionaries representing the query results.
        """
        for feature_name in metadata_feature_names:
            if feature_name not in AVAILABLE_FIELDS_DICT:
                raise exceptions.CellMetadataColumnDoesNotExist(
                    f"Feature {feature_name} is not available for querying. "
                    + f"Please specify any of the following: {', '.join(AVAILABLE_FIELDS_DICT)}."
                )

        if "cas_cell_index" not in metadata_feature_names:
            metadata_feature_names.append("cas_cell_index")

        try:
            return self.cell_operations_dm.get_cell_metadata_by_ids(
                cell_ids=cell_ids, metadata_feature_names=metadata_feature_names
            )
        except dm_exc.NotFound as e:
            raise exceptions.InvalidInputError(str(e))
