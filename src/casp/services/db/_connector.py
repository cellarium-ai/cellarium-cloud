from google.cloud.sql.connector import Connector, IPTypes
import pg8000

import sqlalchemy
from casp import settings


def _connect_with_google_cloud_connector() -> sqlalchemy.engine.base.Engine:
    """
    Initializes a connection pool for a Cloud SQL instance of Postgres.
    Uses the Cloud SQL Python Connector package.

    Refer to GC Documentation:
    https://cloud.google.com/sql/docs/postgres/connect-run#python_1
    """
    ip_type = IPTypes.PRIVATE if settings.DB_PRIVATE_IP else IPTypes.PUBLIC

    # initialize Cloud SQL Python Connector object
    connector = Connector()

    def connection_creator() -> pg8000.dbapi.Connection:
        conn: pg8000.dbapi.Connection = connector.connect(
            instance_connection_string=settings.DB_CONNECTION_NAME,
            driver="pg8000",
            user=settings.DB_USER,
            password=settings.DB_PASSWORD,
            db=settings.DB_NAME,
            ip_type=ip_type,
        )
        return conn

    # The Cloud SQL Python Connector can be used with SQLAlchemy
    # using the 'creator' argument to 'create_engine'
    return sqlalchemy.create_engine(
        "postgresql+pg8000://",
        creator=connection_creator,
        pool_size=5,
        max_overflow=2,
        pool_timeout=30,
        pool_recycle=1800
    )


def _connect_with_regular_db() -> sqlalchemy.engine.base.Engine:
    return sqlalchemy.create_engine(settings.SQLALCHEMY_DATABASE_URI)


def _get_database_engine() -> sqlalchemy.engine.base.Engine:
    if settings.ENVIRONMENT == "local":
        return _connect_with_regular_db()
    elif settings.ENVIRONMENT == "development" or settings.ENVIRONMENT == "production":
        return _connect_with_google_cloud_connector()
    else:
        raise Exception(
            "CAS Database Engine handles one of the following environments: "
            "local, development, production"
        )
