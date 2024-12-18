import io
import json
import os
import random
import tempfile
import typing as t
from datetime import datetime
from unittest.mock import Mock, patch

import anndata
import numpy as np
import pandas as pd
import pytest
from sqlalchemy.orm import sessionmaker

from casp.services import settings
from casp.services.api import services
from casp.services.api.services.cell_operations_service import CellOperationsService
from casp.services.api.services.consensus_engine.strategies.ontology_aware import CellOntologyResource
from casp.services.db import Base, create_engine, models
from casp.services.model_inference.services import ModelInferenceService
from tests.unit import constants, mocks


@pytest.fixture(scope="session", autouse=True)
def patch_google_service_credentials():
    """
    Patch the `get_google_service_credentials` function globally for all tests.
    """
    mock_credentials = Mock()  # Use mock credentials
    mock_project_id = "mock-project-id"

    mock_credentials.token = "mock-token"
    mock_credentials.project_id = mock_project_id
    mock_credentials.get_project_id.return_value = mock_project_id

    def mocked_get_google_service_credentials():
        return mock_credentials, mock_project_id

    # Target the function to be patched
    patch_target = "casp.services.utils.get_google_service_credentials"

    with patch(patch_target, side_effect=mocked_get_google_service_credentials):
        yield


@pytest.fixture(scope="session", autouse=True)
def patch_bigquery_client():
    """
    Patch the BigQuery Client globally for all tests to avoid authentication.
    """
    # Mock the BigQuery client
    mock_client = Mock()

    # Patch the client class to always return the mock instance
    patch_target = "google.cloud.bigquery.Client"
    with patch(patch_target, return_value=mock_client):
        yield mock_client


@pytest.fixture(scope="session", autouse=True)
def test_engine():
    """
    Create a test database engine and apply migrations.
    """
    engine = create_engine()
    Base.metadata.create_all(engine)  # Create tables
    yield engine
    Base.metadata.drop_all(engine)  # Drop tables after all tests


@pytest.fixture(scope="session", autouse=True)
def cleanup_test_db():
    """
    Remove the SQLite test database file after all tests.
    """
    yield  # Let the tests run
    if os.path.exists(settings.TEST_DB_FILE_PATH):
        os.remove(settings.TEST_DB_FILE_PATH)


@pytest.fixture
def cell_info_data() -> t.Dict[int, t.Dict[str, str]]:
    """
    Centralized fixture that represents the source of cell data. It is needed for adding cells into the database and to
    use the same cells for mocking MatchingEngine.

    :return: A dictionary mapping cell IDs (int) to metadata dictionaries with 'cell_type' and
        'cell_type_ontology_term_id'.
    """
    # Define a mapping between cell types and their ontology term IDs
    cell_type_to_ontology_id = {
        "monocyte": "CL:0000576",
        "erythrocyte": "CL:0000232",
        "lymphocyte": "CL:0000542",
        "CD4-positive, alpha-beta T cell": "CL:0000624",
    }

    return {
        i: {
            "cell_type": cell_type,
            "cell_type_ontology_term_id": cell_type_to_ontology_id[cell_type],
        }
        for i, cell_type in enumerate(random.choices(list(cell_type_to_ontology_id.keys()), k=40))
    }


@pytest.fixture
def db_session(test_engine):
    """Provides a session for querying in tests."""
    SessionLocal = sessionmaker(bind=test_engine)
    with SessionLocal() as session:
        yield session


@pytest.fixture(scope="function", autouse=True)
def populate_db(test_engine, db_session, cell_info_data):
    """
    Populate the test database with initial data before each test.
    """
    # Drop all tables and recreate them
    Base.metadata.drop_all(test_engine)
    Base.metadata.create_all(test_engine)
    # Insert sample data
    user_with_quota = models.User(
        id=1,
        username=constants.USER_EMAIL_WITH_QUOTA,
        email=constants.USER_EMAIL_WITH_QUOTA,
        lifetime_cell_quota=constants.USER_LIFETIME_CELL_QUOTA_WITH_QUOTA,
    )
    user_without_quota = models.User(
        id=2,
        username=constants.USER_EMAIL_WITHOUT_QUOTA,
        email=constants.USER_EMAIL_WITHOUT_QUOTA,
        lifetime_cell_quota=constants.USER_LIFETIME_CELL_QUOTA_WITHOUT_QUOTA,
    )
    # Add model
    cas_model = models.CASModel(
        id=1,
        model_name=constants.TEST_MODEL_NAME,
        admin_use_only=constants.TEST_MODEL_ADMIN_USE_ONLY,
        model_file_path=constants.TEST_MODEL_FILE_PATH,
        embedding_dimension=constants.TEST_EMBEDDING_DIMENSION,
    )
    # Add index
    cas_matching_engine_index = models.CASMatchingEngineIndex(
        id=1,
        index_name=constants.TEST_INDEX_NAME,
        model_id=cas_model.id,
        num_neighbors=constants.TEST_INDEX_NUM_NEIGHBORS,
        embedding_dimension=constants.TEST_EMBEDDING_DIMENSION,
        endpoint_id=constants.TEST_INDEX_ENDPOINT_ID,
    )
    # Add User Activity, so that quota of user `user_without_quota` is reached.
    activity1 = models.UserActivity(
        user_id=user_without_quota.id,
        request_id=constants.TEST_REQUEST_ID,
        model_name=constants.TEST_MODEL_NAME,
        method=constants.TEST_USER_ACTIVITY_METHOD,
        cell_count=500,
        event=constants.EVENT_SUCCEEDED_STR,
    )

    activity2 = models.UserActivity(
        user_id=user_without_quota.id,
        request_id=constants.TEST_REQUEST_ID,
        model_name=constants.TEST_MODEL_NAME,
        method=constants.TEST_USER_ACTIVITY_METHOD,
        cell_count=500,
        event=constants.EVENT_SUCCEEDED_STR,
    )

    # Add ingest info into the database
    ingest_info = models.CellIngestInfo(
        cas_ingest_id="ingest_123",
        dataset_id="dataset_001",
        dataset_version_id="v1.0",
        ingest_timestamp=datetime.utcnow(),
    )
    # Add cells and link them to ingest info
    cells = [
        models.CellInfo(
            cas_cell_index=cas_cell_index,
            cell_type=data["cell_type"],  # Extract 'cell_type' from the data
            cell_type_ontology_term_id=data["cell_type_ontology_term_id"],
            cas_ingest_id=ingest_info.cas_ingest_id,  # Link cells to ingest info
        )
        for cas_cell_index, data in cell_info_data.items()
    ]
    db_session.add_all(
        [
            user_with_quota,
            user_without_quota,
            cas_model,
            cas_matching_engine_index,
            activity1,
            activity2,
            ingest_info,
            *cells,
        ]
    )
    db_session.commit()

    # Yield session to allow querying in the test
    yield db_session  # This makes the db_session available for querying during the test


@pytest.fixture
def mock_cell_quota_service() -> services.CellQuotaService:
    """
    Pytest fixture to provide a mocked `CellQuotaDataManager`.

    :return: An instance of `MockCellQuotaDataManager`.
    """
    return services.CellQuotaService(cell_quota_dm=mocks.MockCellQuotaDataManager())


@pytest.fixture
def mock_model_service() -> mocks.MockModelService:
    """
    Fixture for a mocked `ModelInferenceService` with embedding behavior.

    :return: A mocked `ModelInferenceService` instance that generates random embeddings.
    """
    return mocks.MockModelService()


@pytest.fixture
def cell_operations_service_with_mocks(
    mock_cell_quota_service: services.CellQuotaService,
    mock_model_service: ModelInferenceService,
) -> CellOperationsService:
    """
    Fixture to provide a `CellOperationsService` with all mocked dependencies.

    :param mock_cell_quota_service: Mocked `CellQuotaService`.
    :param mock_model_service: Mocked `ModelInferenceService`.

    :return: An instance of `CellOperationsService` with all dependencies mocked.
    """
    service = CellOperationsService(
        # cell_operations_dm=mock_cell_operations_dm,
        # cellarium_general_dm=mock_cellarium_general_dm,
        cell_quota_service=mock_cell_quota_service,
        model_service=mock_model_service,
    )

    return service


@pytest.fixture
def mock_valid_anndata() -> anndata.AnnData:
    """
    Fixture to provide a larger mocked AnnData object.

    :return: An `AnnData` object with random expression data and metadata.
    """
    n_cells = 20  # Number of cells
    n_genes = 500  # Number of genes

    # Generate obs dataframe (metadata for cells)
    obs_data = {
        "total_mrna_umis": np.random.uniform(50.0, 500.0, size=n_cells).astype(np.float32),  # Random float values
    }
    obs = pd.DataFrame(obs_data, index=[f"cell_{i}" for i in range(n_cells)])

    # Generate var dataframe (gene names)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    # Generate X matrix (gene expression data)
    X = np.random.rand(n_cells, n_genes).astype(np.float32)  # Random gene expression values

    # Create AnnData object
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    return adata


@pytest.fixture
def mock_file_with_anndata(mock_valid_anndata: anndata.AnnData) -> io.BytesIO:
    """
    Fixture to serialize the mocked AnnData object into a BytesIO object.

    :param mock_valid_anndata: The AnnData object to serialize.
    :return: A `BytesIO` object containing the serialized AnnData file.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temp_file:
        temp_file_path = temp_file.name
        # Write the AnnData object to a temporary file
        mock_valid_anndata.write_h5ad(temp_file_path)

    # Read the file back into a BytesIO object
    with open(temp_file_path, "rb") as f:
        file = io.BytesIO(f.read())

    return file


@pytest.fixture
def mock_matching_client(cell_info_data: t.Dict[int, t.Any]) -> mocks.MockMatchingClient:
    """
    Fixture for a simplified mocked `MatchingClient` that uses the mock_cell_source.

    :param cell_info_data: Source of cell data with unique IDs and cell types.

    :return: A `MockMatchingClient` instance that generates mock neighbors.
    """
    return mocks.MockMatchingClient(cell_info_data)


def load_ontology_resource_from_file() -> t.Dict[str, t.Any]:
    """
    Load the ontology resource from a predefined JSON file.

    :return: A dictionary representing the simplified ontology resource.
    """
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


@pytest.fixture
def cell_ontology_resource_mock() -> CellOntologyResource:
    """
    Fixture to provide a mocked `CellOntologyResource` instance.

    :return: A `CellOntologyResource` mock instance.
    """
    resource_dict = load_ontology_resource_from_file()
    return CellOntologyResource(cell_ontology_resource_dict=resource_dict)
