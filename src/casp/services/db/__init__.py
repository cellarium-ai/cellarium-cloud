import pg8000
import sqlalchemy.orm
from google.cloud.sql.connector import Connector, IPTypes
from sqlalchemy.orm import declarative_base, sessionmaker

from casp.services import settings

# initialize Cloud SQL Python Connector object
connector = Connector()


def get_private_connection() -> pg8000.dbapi.Connection:
    conn: pg8000.dbapi.Connection = connector.connect(
        instance_connection_string=settings.DB_CONNECTION_NAME,
        driver="pg8000",
        user=settings.DB_USER,
        password=settings.DB_PASSWORD,
        db=settings.DB_NAME,
        enable_iam_auth=False,
        ip_type=IPTypes.PRIVATE,
    )
    return conn


def create_engine() -> sqlalchemy.engine.base.Engine:
    if settings.DB_PRIVATE_IP is not None:
        return sqlalchemy.create_engine(
            url="postgresql+pg8000://",
            creator=get_private_connection,
            pool_size=settings.DB_CONNECTION_POOL_SIZE,
            max_overflow=settings.DB_CONNECTION_POOL_MAX_OVERFLOW,
            pool_timeout=settings.DB_CONNECTION_POOL_TIMEOUT,
            pool_recycle=settings.DB_CONNECTION_POOL_RECYCLE,
        )

    else:
        return sqlalchemy.create_engine(
            url=settings.SQLALCHEMY_DATABASE_URI,
            pool_size=settings.DB_CONNECTION_POOL_SIZE,
            max_overflow=settings.DB_CONNECTION_POOL_MAX_OVERFLOW,
            pool_timeout=settings.DB_CONNECTION_POOL_TIMEOUT,
            pool_recycle=settings.DB_CONNECTION_POOL_RECYCLE,
            echo=settings.ENVIRONMENT == "local",
        )


engine = create_engine()

db_session_maker = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()


def get_db_session_maker() -> sessionmaker:
    import casp.services.db.models  # noqa

    return db_session_maker
