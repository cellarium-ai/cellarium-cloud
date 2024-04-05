"""0010_add_cell_related_models

Revision ID: 7de52712f144
Revises: 310657f0e46c
Create Date: 2024-02-26 16:00:07.051642

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "86eda43c427a"
down_revision = "310657f0e46c"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "cells_cellingestinfo",
        sa.Column("cas_ingest_id", sa.String(length=255), nullable=False),
        sa.Column("dataset_id", sa.String(length=255), nullable=False),
        sa.Column("uns_metadata", sa.JSON(), nullable=False),
        sa.Column("ingest_timestamp", sa.DateTime(), nullable=True),
        sa.PrimaryKeyConstraint("cas_ingest_id"),
    )
    op.create_table(
        "cells_cellfeature",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("cas_feature_index", sa.Integer(), nullable=False),
        sa.Column("original_feature_id", sa.String(length=255), nullable=False),
        sa.Column("feature_name", sa.String(length=255), nullable=False),
        sa.Column("feature_biotype", sa.String(length=255), nullable=True),
        sa.Column("feature_is_filtered", sa.Boolean(), nullable=True),
        sa.Column("feature_reference", sa.String(length=255), nullable=True),
        sa.Column("var_metadata_extra", sa.JSON(), nullable=False),
        sa.Column("cas_ingest_id", sa.String(length=255), nullable=False),
        sa.ForeignKeyConstraint(
            ["cas_ingest_id"],
            ["cells_cellingestinfo.cas_ingest_id"],
        ),
        sa.PrimaryKeyConstraint("cas_feature_index"),
    )
    op.create_table(
        "cells_cellinfo",
        sa.Column("cas_cell_index", sa.Integer(), nullable=False),
        sa.Column("cas_ingest_id", sa.String(length=255), nullable=False),
        sa.Column("original_cell_id", sa.String(length=255), nullable=False),
        sa.Column("obs_metadata_extra", sa.JSON(), nullable=False),
        sa.Column("is_primary_data", sa.Boolean(), nullable=True),
        sa.Column("donor_id", sa.String(length=255), nullable=True),
        sa.Column("cell_type", sa.String(length=255), nullable=False),
        sa.Column("assay", sa.String(length=255), nullable=True),
        sa.Column("development_stage", sa.String(length=255), nullable=True),
        sa.Column("tissue", sa.String(length=255), nullable=True),
        sa.Column("disease", sa.String(length=255), nullable=True),
        sa.Column("organism", sa.String(length=255), nullable=True),
        sa.Column("self_reported_ethnicity", sa.String(length=255), nullable=True),
        sa.Column("sex", sa.String(length=255), nullable=True),
        sa.Column("suspension_type", sa.String(length=255), nullable=True),
        sa.Column("total_mrna_umis", sa.Integer(), nullable=True),
        sa.Column("cell_type_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("assay_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("development_stage_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("tissue_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("disease_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("organism_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("self_reported_ethnicity_ontology_term_id", sa.String(length=255), nullable=True),
        sa.Column("sex_ontology_term_id", sa.String(length=255), nullable=True),
        sa.ForeignKeyConstraint(
            ["cas_ingest_id"],
            ["cells_cellingestinfo.cas_ingest_id"],
        ),
        sa.PrimaryKeyConstraint("cas_cell_index"),
    )


def downgrade() -> None:
    op.drop_table("cells_cellinfo")
    op.drop_table("cells_cellfeature")
    op.drop_table("cells_cellingestinfo")
