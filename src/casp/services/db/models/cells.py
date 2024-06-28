import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import backref, relationship

from casp.services import db


class CellIngestInfo(db.Base):
    cas_ingest_id = sa.Column(sa.String(255), unique=True, nullable=False, primary_key=True)
    dataset_id = sa.Column(sa.String(255), nullable=False)
    dataset_version_id = sa.Column(sa.String(255), nullable=True)
    uns_metadata = sa.Column(JSONB, nullable=False, default={})
    ingest_timestamp = sa.Column(sa.DateTime, default=datetime.datetime.utcnow)

    __tablename__ = "cells_ingestinfo"

    def __str__(self):
        return self.cas_ingest_id


class FeatureInfo(db.Base):
    cas_feature_index = sa.Column(sa.Integer, unique=True, nullable=False, primary_key=True)
    original_feature_id = sa.Column(sa.String(255), unique=False, nullable=False)
    feature_name = sa.Column(sa.String(255), nullable=False)
    feature_biotype = sa.Column(sa.String(255), nullable=True)
    feature_is_filtered = sa.Column(sa.Boolean(), default=False, nullable=True)
    feature_reference = sa.Column(sa.String(255), nullable=True)
    var_metadata_extra = sa.Column(JSONB, nullable=False, default={})
    cas_ingest_id = sa.Column(
        sa.String(255), sa.ForeignKey(f"{CellIngestInfo.__tablename__}.cas_ingest_id"), nullable=False
    )

    __tablename__ = "cells_featureinfo"

    def __str__(self):
        return f"{self.original_feature_id} - {self.feature_name}"


class CellInfo(db.Base):
    # General Information
    cas_cell_index = sa.Column(sa.Integer, unique=True, nullable=False, primary_key=True)
    cas_ingest_id = sa.Column(
        sa.String(255), sa.ForeignKey(f"{CellIngestInfo.__tablename__}.cas_ingest_id"), nullable=False, index=True
    )
    cas_ingest = relationship("CellIngestInfo", backref=backref("cell_info", uselist=False))
    original_cell_id = sa.Column(sa.String(255), nullable=True)
    obs_metadata_extra = sa.Column(JSONB, nullable=True, default={})
    is_primary_data = sa.Column(sa.Boolean(), default=True, nullable=True, index=True)
    donor_id = sa.Column(sa.String(255), nullable=True)
    # Cell Features
    cell_type = sa.Column(sa.String(255), nullable=False, index=True)
    assay = sa.Column(sa.String(255), nullable=True, index=True)
    development_stage = sa.Column(sa.String(255), nullable=True, index=True)
    tissue = sa.Column(sa.String(255), nullable=True, index=True)
    disease = sa.Column(sa.String(255), nullable=True, index=True)
    organism = sa.Column(sa.String(255), nullable=True, index=True)
    self_reported_ethnicity = sa.Column(sa.String(255), nullable=True, index=True)
    sex = sa.Column(sa.String(255), nullable=True, index=True)
    suspension_type = sa.Column(sa.String(255), nullable=True, index=True)
    total_mrna_umis = sa.Column(sa.Integer, nullable=True, index=True)
    # Cell Features Ontology Term IDs
    cell_type_ontology_term_id = sa.Column(sa.String(255), nullable=True, index=True)
    assay_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    development_stage_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    tissue_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    disease_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    organism_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    self_reported_ethnicity_ontology_term_id = sa.Column(sa.String(255), nullable=True)
    sex_ontology_term_id = sa.Column(sa.String(255), nullable=True)

    __tablename__ = "cells_cellinfo"

    def __str__(self):
        return f"{self.cas_cell_index} - {self.original_cell_id}"
