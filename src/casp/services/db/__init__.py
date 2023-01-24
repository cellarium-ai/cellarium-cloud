import sqlalchemy.orm
from sqlalchemy.orm import declarative_base, scoped_session, sessionmaker

from casp import settings

engine = sqlalchemy.create_engine(
    url=settings.SQLALCHEMY_DATABASE_URI,
    pool_size=5,
    max_overflow=2,
    pool_timeout=30,  # 30 seconds
    pool_recycle=1800,  # 30 minutes
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
db_session = scoped_session(SessionLocal)
Base = declarative_base()
Base.query = db_session.query_property()


def save(record) -> None:
    db_session.add(record)
    db_session.commit()


def init_db() -> sqlalchemy.orm.scoped_session:
    import casp.services.db.models  # noqa

    return db_session
