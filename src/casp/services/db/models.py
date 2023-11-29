import datetime

import sqlalchemy as sa
from sqlalchemy.orm import backref, relationship

from casp.services import db, settings


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
    model_name = sa.Column(sa.String(255), unique=True, nullable=False)
    model_file_path = sa.Column(sa.String(255), unique=True, nullable=False)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    admin_use_only = sa.Column(sa.Boolean(), default=True, nullable=False)
    schema_name = sa.Column(sa.String(255), default=settings.DEFAULT_SCHEMA_NAME, nullable=False)
    bq_dataset_name = sa.Column(sa.String(255), default=settings.DEFAULT_MODEL_BQ_DATASET_NAME, nullable=False)
    is_default_model = sa.Column(sa.Boolean(), default=False, nullable=False)
    created_date = sa.Column(sa.DateTime, default=datetime.datetime.utcnow)

    __tablename__ = "cas_model"

    def __str__(self):
        return self.model_name


class CASMatchingEngineIndex(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    index_name = sa.Column(sa.String(255), unique=True, nullable=False)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    endpoint_id = sa.Column(sa.String(255), unique=False, nullable=False)
    deployed_index_id = sa.Column(sa.String(255), unique=True, nullable=False)
    num_neighbors = sa.Column(sa.Integer, nullable=False)
    model_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CASModel.__tablename__}.id"), nullable=False)
    model = relationship("CASModel", backref=backref("cas_matching_engine", uselist=False))

    def __str__(self):
        return self.index_name

    __tablename__ = "cas_matching_engine_index"
