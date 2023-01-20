import sqlalchemy.orm
from sqlalchemy.orm import declarative_base, scoped_session, sessionmaker

from casp.services.db._connector import _get_database_engine

engine = _get_database_engine()
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
