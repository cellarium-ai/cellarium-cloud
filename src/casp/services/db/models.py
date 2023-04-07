import datetime
import sqlalchemy as sa
from sqlalchemy.orm import relationship, backref

from casp.services import db


class User(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    email = sa.Column(sa.String(255), unique=True, nullable=False)
    username = sa.Column(sa.String(255), nullable=False)
    active = sa.Column(sa.Boolean(), default=True, nullable=False)
    requests_processed = sa.Column(sa.Integer, default=0, nullable=False)
    cells_processed = sa.Column(sa.Integer, default=0, nullable=False)
    is_admin = sa.Column(sa.Boolean(), default=True, nullable=False)

    __tablename__ = "user"

    def __repr__(self):
        return self.email


class CASModel(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    system_name = sa.Column(sa.String(255), unique=True, nullable=False)
    model_file_path = sa.Column(sa.String(255), unique=True, nullable=False)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    admin_use_only = sa.Column(sa.Boolean(), default=True, nullable=False)
    created_date = sa.Column(sa.DateTime, default=datetime.datetime.utcnow)

    __tablename__ = "cas_model"


class CASMatchingEngineIndex(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    system_name = sa.Column(sa.String(255), unique=True, nullable=False)
    matching_engine_index_name = sa.Column(sa.String(255), unique=True, nullable=False)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    model_endpoint_url = sa.Column(sa.String(255), unique=True, nullable=True)
    matching_engine_endpoint = sa.Column(sa.String(255), unique=True, nullable=True)
    admin_use_only = sa.Column(sa.Boolean(), default=True, nullable=False)
    cas_model_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CASModel.__tablename__}.id"), nullable=False)
    cas_model = relationship("CASModel", backref=backref("cas_matching_engine_index", uselist=False))

    __tablename__ = "cas_matching_engine_index"
