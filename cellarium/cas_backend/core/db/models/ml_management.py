import datetime

import sqlalchemy as sa
from sqlalchemy.orm import backref, relationship

from cellarium.cas_backend.core import db, settings


class CellInfoMetadata(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    name = sa.Column(sa.String(255), unique=True, nullable=False)
    soma_dataframe_uri = sa.Column(sa.String(1023), nullable=False)
    created_date = sa.Column(sa.DateTime, default=datetime.datetime.utcnow)

    __tablename__ = "ml_management_cell_info_metadata"

    def __str__(self):
        return self.name


class OntologicalColumn(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    cell_info_metadata_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CellInfoMetadata.__tablename__}.id"), nullable=False)
    cell_info_metadata = relationship(
        "CellInfoMetadata", backref=backref("ontological_columns", uselist=True, cascade="all, delete-orphan")
    )
    column_name = sa.Column(sa.String(255), nullable=False)
    ontology_resource_uri = sa.Column(sa.String(1023), nullable=False)
    description = sa.Column(sa.Text, nullable=True)

    __tablename__ = "ml_management_ontological_column"
    __table_args__ = (sa.UniqueConstraint("cell_info_metadata_id", "column_name"),)

    def __str__(self):
        return self.column_name


class CASModel(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    model_name = sa.Column(sa.String(255), unique=True, nullable=False)
    model_file_path = sa.Column(sa.String(255), unique=False, nullable=False)
    description = sa.Column(sa.String(511), nullable=True)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    admin_use_only = sa.Column(sa.Boolean(), default=True, nullable=False)
    schema_name = sa.Column(sa.String(255), default=settings.DEFAULT_SCHEMA_NAME, nullable=False)
    bq_dataset_name = sa.Column(sa.String(255), default=settings.DEFAULT_MODEL_BQ_DATASET_NAME, nullable=False)
    is_default_model = sa.Column(sa.Boolean(), default=False, nullable=False)
    cell_info_metadata_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CellInfoMetadata.__tablename__}.id"), nullable=False)
    cell_info_metadata = relationship("CellInfoMetadata", backref=backref("models", uselist=True))
    created_date = sa.Column(sa.DateTime, default=datetime.datetime.utcnow)

    __tablename__ = "ml_management_model"

    def __str__(self):
        return self.model_name


class CASMatchingEngineIndex(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    index_name = sa.Column(sa.String(255), unique=True, nullable=False)
    description = sa.Column(sa.String(511), nullable=True)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    endpoint_id = sa.Column(sa.String(255), unique=False, nullable=False)
    deployed_index_id = sa.Column(sa.String(255), unique=True, nullable=True)
    num_neighbors = sa.Column(sa.Integer, nullable=False)
    model_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CASModel.__tablename__}.id"), nullable=False)
    model = relationship("CASModel", backref=backref("cas_matching_engine", uselist=False))
    is_grpc = sa.Column(sa.Boolean(), default=True, nullable=False)
    api_endpoint = sa.Column(sa.String(255), nullable=True)

    def __str__(self):
        return self.index_name

    __tablename__ = "ml_management_matchingengineindex"


class CASVectorIndex(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    index_name = sa.Column(sa.String(255), unique=True, nullable=False)
    description = sa.Column(sa.String(511), nullable=True)
    embedding_dimension = sa.Column(sa.Integer, nullable=False)
    num_neighbors = sa.Column(sa.Integer, nullable=False)
    index_uri = sa.Column(sa.String(1023), nullable=False)
    index_type = sa.Column(sa.String(255), nullable=False)
    distance_metric = sa.Column(sa.String(255), nullable=False)
    nprobe = sa.Column(sa.Integer, nullable=True)
    l_search = sa.Column(sa.Integer, nullable=True)
    model_id = sa.Column(sa.Integer, sa.ForeignKey(f"{CASModel.__tablename__}.id"), nullable=False, unique=True)
    model = relationship("CASModel", backref=backref("cas_vector_index", uselist=False))

    def __str__(self):
        return self.index_name

    __tablename__ = "ml_management_vectorindex"
