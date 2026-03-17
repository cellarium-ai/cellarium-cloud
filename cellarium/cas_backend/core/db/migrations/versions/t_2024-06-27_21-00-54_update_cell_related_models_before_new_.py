"""Update cell related models before new ingest

Revision ID: 51a4ad2bb22d
Revises: 05b03bfc2e69
Create Date: 2024-06-27 21:00:54.062701

"""

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "51a4ad2bb22d"
down_revision = "05b03bfc2e69"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.alter_column("cells_cellinfo", "original_cell_id", existing_type=sa.VARCHAR(length=255), nullable=True)
    op.alter_column(
        "cells_cellinfo", "obs_metadata_extra", existing_type=postgresql.JSONB(astext_type=sa.Text()), nullable=True
    )
    op.add_column("cells_ingestinfo", sa.Column("dataset_version_id", sa.String(length=255), nullable=True))


def downgrade() -> None:
    op.drop_column("cells_ingestinfo", "dataset_version_id")
    op.alter_column(
        "cells_cellinfo", "obs_metadata_extra", existing_type=postgresql.JSONB(astext_type=sa.Text()), nullable=False
    )
    op.alter_column("cells_cellinfo", "original_cell_id", existing_type=sa.VARCHAR(length=255), nullable=False)
