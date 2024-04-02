import sqlalchemy.orm
from sqlalchemy.orm import declarative_base, scoped_session, sessionmaker

from casp.services import settings

engine = sqlalchemy.create_engine(
    url=settings.SQLALCHEMY_DATABASE_URI,
    pool_size=settings.DB_CONNECTION_POOL_SIZE,
    max_overflow=settings.DB_CONNECTION_POOL_MAX_OVERFLOW,
    pool_timeout=settings.DB_CONNECTION_POOL_TIMEOUT,
    pool_recycle=settings.DB_CONNECTION_POOL_RECYCLE,
)
db_session_maker = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()



def get_db_session_maker() -> sessionmaker:
    import casp.services.db.models  # noqa

    return db_session_maker
