"""
Fixtures for setting up and managing the database and test environment for unit tests.

This module provides a set of pytest fixtures for initializing and interacting with the test database,
resetting the database schema between tests, and populating the database with initial data for testing.
The fixtures also include necessary session management and integration with SQLAlchemy.
"""

import os
import typing as t
from datetime import datetime

import pytest
import sqlalchemy as sa
from sqlalchemy.orm import Session, sessionmaker

from casp.services import settings
from casp.services.api import schemas
from casp.services.db import Base, create_engine, models
from tests.unit.fixtures import constants


@pytest.fixture(scope="session", autouse=True)
def test_engine() -> t.Generator[sa.engine.Engine, None, None]:
    """
    Fixture to create a test database engine. Create all tables before tests and remove tables and test database file
    after

    **Scope**: Session

    **Autouse**: Yes

    :return: A SQLAlchemy engine connected to the test database.
    """
    # Execute before all tests
    engine = create_engine()
    Base.metadata.create_all(engine)
    yield engine  # Let the tests run
    # Execute after all tests
    Base.metadata.drop_all(engine)

    if os.path.exists(settings.TEST_DB_FILE_PATH):
        os.remove(settings.TEST_DB_FILE_PATH)


@pytest.fixture
def db_session(test_engine: sa.engine.Engine) -> t.Generator[Session, None, None]:
    """
    Fixture to provide a SQLAlchemy session for querying during tests.

    **Scope**: Function

    :param test_engine: The SQLAlchemy engine used for the test database.

    :return: A generator that yields a SQLAlchemy session.
    """
    session_local_class = sessionmaker(bind=test_engine)
    with session_local_class() as session:
        yield session


@pytest.fixture
def reset_db(test_engine: sa.engine.Engine) -> t.Generator[None, None, None]:
    """
    Fixture to reset the database schema by dropping and recreating all tables.

    **Scope**: Function

    **Autouse**: False

    :param test_engine: The SQLAlchemy engine used for the test database.
    """
    Base.metadata.drop_all(test_engine)
    Base.metadata.create_all(test_engine)
    yield


@pytest.fixture
def populate_db(
    db_session: Session, reset_db: None, cell_info_data: t.List[schemas.CellariumCellMetadata]
) -> t.Generator[None, None, None]:
    """
    Fixture to populate the database with initial data for tests. The database schema is reset before adding the data.

    **Scope**: Function

    **Autouse**: False

    :param db_session: The SQLAlchemy session used for database operations.
    :param reset_db: Fixture to reset the database schema.
    :param cell_info_data: List of cell metadata to be inserted into the database.
    """
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
    activity_1 = models.UserActivity(
        user_id=user_without_quota.id,
        request_id=constants.TEST_REQUEST_ID,
        model_name=constants.TEST_MODEL_NAME,
        method=constants.TEST_USER_ACTIVITY_METHOD,
        cell_count=constants.USER_LIFETIME_CELL_QUOTA_WITHOUT_QUOTA,
        event=constants.EVENT_SUCCEEDED_STR,
    )

    activity_2 = models.UserActivity(
        user_id=user_without_quota.id,
        request_id=constants.TEST_REQUEST_ID,
        model_name=constants.TEST_MODEL_NAME,
        method=constants.TEST_USER_ACTIVITY_METHOD,
        cell_count=constants.USER_LIFETIME_CELL_QUOTA_WITHOUT_QUOTA,
        event=constants.EVENT_SUCCEEDED_STR,
    )

    # Add ingest info into the database
    ingest_info = models.CellIngestInfo(
        cas_ingest_id=constants.TEST_CAS_INGEST_ID,
        dataset_id=constants.TEST_CAS_DATASET_ID,
        dataset_version_id=constants.TEST_CAS_DATASET_VERSION_ID,
        ingest_timestamp=datetime.utcnow(),
    )
    # Add cells and link them to ingest info
    cells = [
        models.CellInfo(
            cas_cell_index=cell_info.cas_cell_index,
            cell_type=cell_info.cell_type,
            cell_type_ontology_term_id=cell_info.cell_type_ontology_term_id,
            cas_ingest_id=ingest_info.cas_ingest_id,
        )
        for cell_info in cell_info_data
    ]
    db_session.add_all(
        [
            user_with_quota,
            user_without_quota,
            cas_model,
            cas_matching_engine_index,
            activity_1,
            activity_2,
            ingest_info,
            *cells,
        ]
    )
    db_session.commit()

    yield
