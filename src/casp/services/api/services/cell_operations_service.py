"""
Cellarium Service Controller. It provides methods to communicate with services in Cellarium Cloud
infrastructure over different protocols in async manner.
"""

import logging
import math
import typing as t

import anndata
import numpy as np
from tenacity import AsyncRetrying, stop_after_attempt, wait_exponential, wait_random

from casp.services import settings
from casp.services.api import clients, schemas
from casp.services.api.clients.matching_client import MatchingClient, MatchResult
from casp.services.api.data_manager import CellariumGeneralDataManager, CellOperationsDataManager
from casp.services.api.data_manager import exceptions as dm_exc
from casp.services.api.services import consensus_engine, exceptions
from casp.services.db import models
from casp.services.utils import numpy_utils

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
        cell_operations_dm: CellOperationsDataManager = None,
        cellarium_general_dm: CellariumGeneralDataManager = None,
    ):
        self.cell_operations_dm = cell_operations_dm or CellOperationsDataManager()
        self.cellarium_general_dm = cellarium_general_dm or CellariumGeneralDataManager()

    def authorize_model_for_user(self, user: models.User, model_name: str) -> models.CASModel:
        """
        Authorize user to use a specific model. If user is not authorized, raise an exception.

        :param user: User object to check permissions for.
        :param model_name: Model name to check permissions for.

        :return: Model object if user is authorized.
        """
        try:
            model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)
        except dm_exc.NotFound as e:
            raise exceptions.InvalidInputError(str(e))
        if model.admin_use_only and not user.is_admin:
            raise exceptions.AccessDeniedError(
                f"{model_name} model is not available. Please reach out to the Cellarium team for more information."
            )
        return model

    def __get_cell_count_from_anndata(self, file: t.BinaryIO) -> int:
        """
        Get the number of cells in the anndata file.

        :param file: Anndata file to get the number of cells from.

        :return: Number of cells in the anndata file.
        """
        adata = anndata.read_h5ad(file)
        cell_count = len(adata)
        # Reset the file pointer to the beginning of the file because we're gonna need to read it
        # again later.
        file.seek(0)
        return cell_count

    def __verify_quota_and_log_activity(self, user: models.User, file: t.BinaryIO, model_name: str, method: str) -> int:
        """
        Verify that the number of cells in the anndata file does not exceed the user's remaining
        quota and log the user activity in the database in an atomic transaction.

        :param user: User object to check quota for.
        :param file: Anndata file to check the number of cells in.
        :param model_name: Model name to use for annotation.
        :param method: Method name to log in the database.

        :return: Number of cells in the anndata file.

        :raises exceptions.QuotaExceededException: If the number of cells in the anndata file exceeds the user's quota.
        """
        cell_count = self.__get_cell_count_from_anndata(file=file)
        user_quota = self.cellarium_general_dm.get_remaining_quota_for_user(user=user)
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

    def __get_match_index_endpoint_for_model(self, model_name: str) -> models.CASMatchingEngineIndex:
        return self.cellarium_general_dm.get_index_for_model(model_name=model_name)

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
        self, embeddings: np.array, model_name: str, chunk_matches_function: t.Callable
    ) -> MatchResult:
        """
        Get KNN matches for embeddings broken into chunks with retry logic

        :param embeddings: Embeddings to match.
        :param model_name: Model name to use for matching.

        :return: List of lists of MatchNeighbor objects.
        """
        index = self.__get_match_index_endpoint_for_model(model_name=model_name)

        matching_client = MatchingClient.from_index(index)

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

    async def get_knn_matches_from_embeddings(self, embeddings: np.array, model_name: str) -> MatchResult:
        """
        Run KNN matching synchronously using Matching Engine client over gRPC.

        :param embeddings: Embeddings from the model.
        :param model_name: Model name to use for matching.

        :return: Matches in a MatchResult object
        """

        return await self.__get_knn_matches_from_embeddings_with_retry(
            embeddings=embeddings,
            model_name=model_name,
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
        num_chunks: int = math.ceil(len(embeddings) / chunk_size)
        return np.array_split(embeddings, num_chunks)

    async def get_knn_matches(self, file: t.Any, model_name: str) -> t.Tuple[t.List[str], MatchResult]:
        """
        Get KNN matches for anndata object.

        :param file: Anndata object to get KNN matches for.
        :param model_name: Model name to use for matching.

        :return: Matches in a MatchResult object
        """
        logger.info("Getting embeddings")
        query_ids, embeddings = await CellOperationsService.get_embeddings(file_to_embed=file, model_name=model_name)

        if embeddings.size == 0:
            # No further processing needed if there are no embeddings, return empty match result
            return [], MatchResult()

        logger.info("Getting KNN matches")
        knn_response = await self.get_knn_matches_from_embeddings(embeddings=embeddings, model_name=model_name)

        return query_ids, knn_response

    async def annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        self,
        user: models.User,
        file: t.BinaryIO,
        model_name: str,
        include_extended_output: t.Optional[bool] = None,
    ) -> schemas.QueryAnnotationCellTypeSummaryStatisticsType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.  This is a wrapper method
        for annotate_cell_type_summary_statistics_strategy that logs to the user activity table and verifies the input
        file doesn't exceed the user's remaining quota

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
        :param include_extended_output: Boolean flag indicating whether to include dev metadata in the response. Used only
            in `cell_type_count` method.

        :return: JSON response with annotations.
        """
        # Make sure the user has enough quota to process the file
        cell_count = self.__verify_quota_and_log_activity(
            user=user, file=file, model_name=model_name, method="annotate_cell_type_summary_statistics_strategy"
        )

        # Annotate the file and log successful activity (or log failure if an exception is raised)
        try:
            annotation_response = await self.annotate_cell_type_summary_statistics_strategy(
                user=user, file=file, model_name=model_name, include_extended_output=include_extended_output
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

    async def annotate_cell_type_summary_statistics_strategy(
        self,
        user: models.User,
        file: t.BinaryIO,
        model_name: str,
        include_extended_output: t.Optional[bool] = None,
    ) -> schemas.QueryAnnotationCellTypeSummaryStatisticsType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
        :param include_extended_output: Boolean flag indicating whether to include dev metadata in the response. Used only
            in `cell_type_count` method.

        :return: JSON response with annotations.
        """
        cas_model = self.authorize_model_for_user(user=user, model_name=model_name)
        query_ids, knn_response = await self.get_knn_matches(file=file, model_name=model_name)

        strategy = consensus_engine.CellTypeSummaryStatisticsConsensusStrategy(
            cell_operations_dm=self.cell_operations_dm,
            cas_model=cas_model,
            include_extended_output=include_extended_output,
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
        # Make sure the user has enough quota to process the file
        cell_count = self.__verify_quota_and_log_activity(
            user=user, file=file, model_name=model_name, method="annotate_cell_type_ontology_aware_strategy"
        )

        # Annotate the file and log successful activity (or log failure if an exception is raised)
        try:
            annotation_response = await self.annotate_cell_type_ontology_aware_strategy(
                user=user,
                file=file,
                model_name=model_name,
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

    async def annotate_cell_type_ontology_aware_strategy(
        self, user: models.User, file: t.BinaryIO, model_name: str, prune_threshold: float, weighting_prefactor: float
    ) -> schemas.QueryAnnotationOntologyAwareType:
        """
        Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
        model schema. Increment user cells processed counter after successful annotation.

        :param user: User object used to increment user cells processed counter.
        :param file: Byte object of :class:`anndata.AnnData` file to annotate.
        :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
        :param prune_threshold: Prune threshold for the ontology-aware annotation strategy.
        :param weighting_prefactor: Distance exponential weighting prefactor.

        :return: JSON response with annotations.
        """
        _ = self.authorize_model_for_user(user=user, model_name=model_name)
        query_ids, knn_response = await self.get_knn_matches(file=file, model_name=model_name)

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
        self.authorize_model_for_user(user=user, model_name=model_name)
        query_ids, embeddings = await CellOperationsService.get_embeddings(file_to_embed=file, model_name=model_name)
        if embeddings.size == 0:
            # No further processing needed if there are no embeddings
            return []
        knn_response = await self.get_knn_matches_from_embeddings(embeddings=embeddings, model_name=model_name)

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
        self, user: models.User, cell_ids: t.List[int], metadata_feature_names: t.List[str], model_name: str
    ) -> t.List[schemas.CellariumCellMetadata]:
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
