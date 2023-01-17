import sqlalchemy.orm
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker, scoped_session
from casp.admin import settings

engine = create_engine(settings.SQLALCHEMY_DATABASE_URI)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
db_session = scoped_session(SessionLocal)
Base = declarative_base()
Base.query = db_session.query_property()


def save(record) -> None:
    db_session.add(record)
    db_session.commit()


def init_db() -> sqlalchemy.orm.scoped_session:
    import casp.db.models
    return db_session


def create_all():
    Base.metadata.create_all(bind=engine)
