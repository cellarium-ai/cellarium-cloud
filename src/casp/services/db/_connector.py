from google.cloud.sql.connector import Connector, IPTypes
import pg8000

import sqlalchemy
from casp import settings


def connect_unix_socket() -> sqlalchemy.engine.base.Engine:
    """ Initializes a Unix socket connection pool for a Cloud SQL instance of Postgres. """
    # Note: Saving credentials in environment variables is convenient, but not
    # secure - consider a more secure solution such as
    # Cloud Secret Manager (https://cloud.google.com/secret-manager) to help
    # keep secrets safe.
    db_user = settings.DB_USER
    db_pass = settings.DB_PASSWORD
    db_name = settings.DB_NAME
    unix_socket_path = settings.DB_INSTANCE_UNIX_SOCKET
    print("I AM IN UNIX SOCKET INIT")
    pool = sqlalchemy.create_engine(
        sqlalchemy.engine.url.URL.create(
            drivername="postgresql+pg8000",
            username=db_user,
            password=db_pass,
            database=db_name,
            query={"unix_sock": "{}/.s.PGSQL.5432".format(unix_socket_path)},
        ),
        pool_size=5,
        max_overflow=2,
        pool_timeout=30,  # 30 seconds
        pool_recycle=1800,  # 30 minutes
    )
    return pool


def _init_regular_connection_engine() -> sqlalchemy.engine.base.Engine:
    print("I AM IN INIT REGULAR CONNECTION")
    return sqlalchemy.create_engine(settings.SQLALCHEMY_DATABASE_URI)


def _get_database_engine() -> sqlalchemy.engine.base.Engine:
    return connect_unix_socket()
    # if settings.ENVIRONMENT == "local":
    #     return _init_regular_connection_engine()
    # elif settings.ENVIRONMENT == "development" or settings.ENVIRONMENT == "production":
    #     return connect_unix_socket()
    # else:
    #     raise Exception(
    #         "CAS Database Engine handles one of the following environments: "
    #         "local, development, production"
    #     )
