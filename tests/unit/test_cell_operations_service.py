"""
Tests things in the cell_operations_service that can reasonably be tested with a unit test.
"""

import io
import re
import typing as t

import anndata
import numpy
import numpy as np
import pytest
from mockito import matchers, mock, unstub, verify, when
from parameterized import parameterized

from casp.services.api.clients.matching_client import MatchingClient, MatchResult
from casp.services.api.services import exceptions
from casp.services.api.services.cell_operations_service import CellOperationsService
from casp.services.db import models
from tests.unit.test_utils import async_return

USER_ADMIN = models.User(id=1, is_admin=True)
USER_NON_ADMIN = models.User(id=2, is_admin=False)
MODEL = models.CASModel(id=1, model_name="model_name", admin_use_only=False)
MODEL_ADMIN_ONLY = models.CASModel(id=2, model_name="admin_only_model", admin_use_only=True)
INDEX = models.CASMatchingEngineIndex(
    is_grpc=True,
    endpoint_id="endpoint_id_grpc",
    deployed_index_id="deployed_index_id_grpc",
    num_neighbors=3,
    model_id=MODEL.id,
    model=MODEL,
)
ANNDATA_DATA = b"testdata"
ANNDATA_FILE = io.BytesIO(ANNDATA_DATA)


class TestCellOperationsService:
    """
    Test the CellOperationsService class.

    """

    def setup_method(self) -> None:
        self.cell_operations_service = CellOperationsService(
            cell_operations_dm=mock(), cellarium_general_dm=mock(), cell_quota_service=mock()
        )

    def teardown_method(self) -> None:
        unstub()

    def test_split_embeddings_into_chunks_even_split(self) -> None:
        embeddings = numpy.random.rand(100, 100)
        chunk_size = 10

        # Turns out this is how you call a private method in a Python unit test.
        chunks = CellOperationsService._CellOperationsService__split_embeddings_into_chunks(embeddings, chunk_size)

        assert len(chunks) == 10
        assert sum(len(chunk) for chunk in chunks) == 100

    def test_split_embeddings_into_chunks_uneven_split(self) -> None:
        embeddings = numpy.random.rand(100, 100)
        chunk_size = 9

        # Turns out this is how you call a private method in a Python unit test.
        chunks = CellOperationsService._CellOperationsService__split_embeddings_into_chunks(embeddings, chunk_size)

        assert len(chunks) == 12
        assert sum(len(chunk) for chunk in chunks) == 100

    def test_mismatched_embeddings_and_queries(self) -> None:
        embeddings = [[0, 1, 2]]
        knn_response = MatchResult(matches=[])
        self.__mock_apis(
            model=MODEL,
            index=INDEX,
            adata=ANNDATA_DATA,
            embeddings=embeddings,
            matching_client_response=knn_response,
        )
        with pytest.raises(
            exceptions.VectorSearchResponseError,
            match=re.escape("Number of query ids (1) and knn matches (0) does not match."),
        ):
            CellOperationsService._CellOperationsService__validate_knn_response(
                embeddings=embeddings, knn_response=knn_response
            )

    def test_empty_neighbors(self) -> None:
        embeddings = [[0, 1, 2]]
        knn_response = MatchResult(matches=[MatchResult.NearestNeighbors(neighbors=[])])
        self.__mock_apis(
            model=MODEL,
            index=INDEX,
            adata=ANNDATA_DATA,
            embeddings=embeddings,
            matching_client_response=knn_response,
        )
        with pytest.raises(
            exceptions.VectorSearchResponseError, match=re.escape("Vector Search returned a match with 0 neighbors.")
        ):
            CellOperationsService._CellOperationsService__validate_knn_response(
                embeddings=embeddings, knn_response=knn_response
            )

    @parameterized.expand(
        [
            ([], [], False),
            ([[0, 1, 2]], ["erythrocyte"], False),
            ([[0, 1, 2], [3, 4, 5]], ["erythrocyte", "monocyte"], True),
        ]
    )
    @pytest.mark.asyncio
    async def test_annotate_adata_file(
        self, embeddings: t.List[t.List[float]], cell_types: t.List[str], include_dev_metadata: bool
    ) -> None:
        """
        Test the annotate_adata_file method.

        :param embeddings: The embeddings to use for the test.
        :param cell_types: The cell types to use for the test.  Should match 1:1 with the embeddings array.
        :param include_dev_metadata: Whether or not to set the include_dev_metadata flag when calling the method.
        """
        matching_client_response = self.__mock_apis(model=MODEL, index=INDEX, adata=embeddings, embeddings=embeddings)

        # mock calls to get cell distribution
        temp_table_fqn = "temp_table_fqn"
        query_ids = [f"q{i}" for i in range(len(embeddings))]
        response = [
            {"query_cell_id": query_ids[i], "matches": [{"cell_type": cell_types[i], "cell_count": 10}]}
            for i in range(len(embeddings))
        ]
        when(self.cell_operations_service.cell_operations_dm).insert_matches_to_temp_table(
            query_ids=query_ids, knn_response=matching_client_response
        ).thenReturn(temp_table_fqn)
        if include_dev_metadata:
            when(self.cell_operations_service.cell_operations_dm).get_neighborhood_distance_summary_dev_details(
                cas_model=MODEL, match_temp_table_fqn=temp_table_fqn
            ).thenReturn(response)
        else:
            when(self.cell_operations_service.cell_operations_dm).get_neighborhood_distance_summary(
                cas_model=MODEL, match_temp_table_fqn=temp_table_fqn
            ).thenReturn(response)

        when(self.cell_operations_service)._CellOperationsService__read_and_validate_anndata_file(
            file=ANNDATA_FILE
        ).thenReturn(embeddings)
        when(self.cell_operations_service.cell_quota_service).get_remaining_quota_for_user(user=USER_ADMIN).thenReturn(
            10000
        )

        actual_response = (
            await self.cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
                user=USER_ADMIN,
                file=ANNDATA_FILE,
                model_name=MODEL.model_name,
                include_extended_output=include_dev_metadata,
            )
        )
        assert actual_response == response

        verify(self.cell_operations_service.cellarium_general_dm).log_user_activity(
            user_id=USER_ADMIN.id,
            model_name=MODEL.model_name,
            method="annotate_cell_type_summary_statistics_strategy",
            cell_count=len(query_ids),
            event=models.UserActivityEvent.STARTED,
        )

        verify(self.cell_operations_service.cellarium_general_dm).log_user_activity(
            user_id=USER_ADMIN.id,
            model_name=MODEL.model_name,
            method="annotate_cell_type_summary_statistics_strategy",
            cell_count=len(query_ids),
            event=models.UserActivityEvent.SUCCEEDED,
        )

    @pytest.mark.asyncio
    async def test_annotate_adata_file_quota_exceeded(self) -> None:
        """
        Test the annotate_adata_file method returns the correct error in the case of an exceeded
        quota.
        """
        self.__mock_apis(model=MODEL, index=INDEX, adata=ANNDATA_DATA, embeddings=[])

        when(self.cell_operations_service)._CellOperationsService__read_and_validate_anndata_file(
            file=ANNDATA_FILE
        ).thenReturn(ANNDATA_DATA)
        when(self.cell_operations_service.cell_quota_service).get_remaining_quota_for_user(user=USER_ADMIN).thenReturn(
            2
        )

        with pytest.raises(exceptions.QuotaExceededException):
            await self.cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
                user=USER_ADMIN,
                file=ANNDATA_FILE,
                model_name=MODEL.model_name,
                include_extended_output=False,
            )

    @pytest.mark.asyncio
    async def test_annotate_adata_file_invalid_anndata_dtype(self) -> None:
        """
        Test the annotate_adata_file method returns the correct error in the case of the anndata
        matrix being of the wrong dtype
        """
        self.__mock_apis(model=MODEL, index=INDEX, adata=ANNDATA_DATA, embeddings=[])

        mock_anndata = anndata.AnnData(np.array([[0, 1, 2]]))
        mock_anndata.X = mock_anndata.X.astype(np.float64)

        when(anndata).read_h5ad(ANNDATA_FILE).thenReturn(mock_anndata)

        with pytest.raises(exceptions.InvalidAnndataFileException):
            await self.cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
                user=USER_ADMIN,
                file=ANNDATA_FILE,
                model_name=MODEL.model_name,
                include_extended_output=False,
            )

    @parameterized.expand(
        [
            ([], []),
            ([[0, 1, 2]], [{"query_cell_id": "q0", "neighbors": [{"cas_cell_index": "0", "distance": 0.0}]}]),
            (
                [[0, 1, 2], [3, 4, 5]],
                [
                    {"query_cell_id": "q0", "neighbors": [{"cas_cell_index": "0", "distance": 0.0}]},
                    {"query_cell_id": "q1", "neighbors": [{"cas_cell_index": "1", "distance": 0.0}]},
                ],
            ),
        ]
    )
    @pytest.mark.asyncio
    async def test_search_adata_file(
        self, embeddings: t.List[t.List[float]], expected_response: t.List[t.Dict[str, t.Any]]
    ) -> None:
        """
        Test the search_adata_file method.

        :param embeddings: The embeddings to use for the test.
        :param expected_response: The expected response from the method.

        """
        self.__mock_apis(model=MODEL, index=INDEX, adata=ANNDATA_DATA, embeddings=embeddings)

        when(self.cell_operations_service)._CellOperationsService__read_and_validate_anndata_file(
            file=ANNDATA_FILE
        ).thenReturn(ANNDATA_DATA)

        actual_response = await self.cell_operations_service.search_adata_file(
            user=USER_ADMIN, file=ANNDATA_FILE, model_name=MODEL.model_name
        )
        assert actual_response == expected_response

    @pytest.mark.asyncio
    async def test_get_cells_by_ids_for_user(self) -> None:
        self.__mock_apis()
        cell_ids = [1, 2, 3]
        metadata_feature_names = ["cell_type", "assay"]

        self.cell_operations_service.get_cells_by_ids_for_user(
            cell_ids=cell_ids,
            metadata_feature_names=metadata_feature_names,
        )
        verify(self.cell_operations_service.cell_operations_dm).get_cell_metadata_by_ids(
            cell_ids=cell_ids,
            metadata_feature_names=metadata_feature_names,
        )

    @pytest.mark.asyncio
    async def test_get_cells_by_ids_for_user_bad_feature_name(self) -> None:
        self.__mock_apis()
        cell_ids = [1, 2, 3]
        metadata_feature_names = ["cell_type", "foo"]

        with pytest.raises(
            exceptions.CellMetadataColumnDoesNotExist, match=re.escape("Feature foo is not available for querying.")
        ):
            self.cell_operations_service.get_cells_by_ids_for_user(
                cell_ids=cell_ids,
                metadata_feature_names=metadata_feature_names,
            )

    def __mock_apis(
        self,
        model: models.CASModel = MODEL,
        index: models.CASMatchingEngineIndex = INDEX,
        adata: anndata.AnnData = ANNDATA_DATA,
        embeddings: t.List[t.List[float]] = [],
        matching_client_response: t.Optional[MatchResult] = None,
    ) -> MatchResult:
        """
        Mock call to the model embedding service and the matching client.

        :param model: The model to mock.
        :param index: The index to mock.  This should be an index for the model.
        :param anndata_data: The source anndata to mock
        :param embeddings: The returned embeddings to mock.
        :param matching_client_response: The response from the matching client. If it isn't passed in,
        a response will be created based on the embeddings.

        :return: The response from the matching client to be used for further mocking.
        """

        # mock calls to model embedding service
        model_name = model.model_name
        embeddings = np.array(embeddings, dtype=np.float32)
        query_ids = [f"q{i}" for i in range(len(embeddings))]
        when(self.cell_operations_service.model_service).embed_adata(adata=adata, model=model).thenReturn(
            (query_ids, embeddings)
        )

        when(self.cell_operations_service.cellarium_general_dm).get_model_by_name(model_name=model_name).thenReturn(
            model
        )

        # mock calls to the matching client
        matching_client = mock()
        matching_client_response = (
            MatchResult(
                matches=[
                    MatchResult.NearestNeighbors(
                        neighbors=[
                            MatchResult.Neighbor(cas_cell_index=str(i), distance=0.0, feature_vector=list(embedding))
                        ]
                    )
                    for i, embedding in enumerate(embeddings)
                ]
            )
            if matching_client_response is None
            else matching_client_response
        )

        # Note: comparing numpy arrays with == doesn't work well with mockito so we need to be a little clever here
        when(matching_client).match(queries=matchers.arg_that(lambda arg: np.array_equal(arg, embeddings))).thenReturn(
            async_return(matching_client_response)
        )
        when(MatchingClient).from_index(index).thenReturn(matching_client)

        return matching_client_response
