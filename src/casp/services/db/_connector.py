from google.cloud.sql.connector import Connector, IPTypes
import pg8000

import sqlalchemy
from casp import settings


# def _connect_with_google_cloud_connector() -> sqlalchemy.engine.base.Engine:
#     """
#     Initializes a connection pool for a Cloud SQL instance of Postgres.
#     Uses the Cloud SQL Python Connector package.
#
#     Refer to GC Documentation:
#     https://cloud.google.com/sql/docs/postgres/connect-run#python_1
#     """
#
#     # initialize Cloud SQL Python Connector object
#     connector = Connector()
#
#     def connection_creator() -> pg8000.dbapi.Connection:
#         conn: pg8000.dbapi.Connection = connector.connect(
#             instance_connection_string=settings.DB_CONNECTION_NAME,
#             driver="pg8000",
#             user=settings.DB_USER,
#             password=settings.DB_PASSWORD,
#             db=settings.DB_NAME,
#             ip_type=IPTypes.PUBLIC,
#         )
#         return conn
#
#     # The Cloud SQL Python Connector can be used with SQLAlchemy
#     # using the 'creator' argument to 'create_engine'
#     return sqlalchemy.create_engine(
#         "postgresql+pg8000://",
#         creator=connection_creator,
#         pool_size=5,
#         max_overflow=2,
#         pool_timeout=30,
#         pool_recycle=1800
#     )


def _init_db_connection_unix():
    db_config = {
        'pool_size': 5,
        'max_overflow': 2,
        'pool_timeout': 30,
        'pool_recycle': 1800,
    }
    return _init_unix_connection_engine(db_config)


def _init_unix_connection_engine(db_config):
    pool = sqlalchemy.create_engine(
        sqlalchemy.engine.url.URL(
            drivername="postgres+pg8000",
            host=settings.DB_HOST,
            port=settings.DB_PORT,
            username=settings.DB_USER,
            password=settings.DB_PASSWORD,
            database=settings.DB_NAME,
        ),
        **db_config
    )
    pool.dialect.description_encoding = None
    return pool


def _init_regular_connection_engine() -> sqlalchemy.engine.base.Engine:
    return sqlalchemy.create_engine(settings.SQLALCHEMY_DATABASE_URI)


def _get_database_engine() -> sqlalchemy.engine.base.Engine:
    if settings.ENVIRONMENT == "local":
        return _init_regular_connection_engine()
    elif settings.ENVIRONMENT == "development" or settings.ENVIRONMENT == "production":
        return _init_db_connection_unix()
    else:
        raise Exception(
            "CAS Database Engine handles one of the following environments: "
            "local, development, production"
        )
