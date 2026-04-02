"""add_cell_metadata_uri_to_model

Revision ID: a1b2c3d4e5f6
Revises: 532fd37ec55a
Create Date: 2026-03-30 00:00:00.000000

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "a1b2c3d4e5f6"
down_revision = "532fd37ec55a"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "ml_management_model",
        sa.Column("cell_metadata_uri", sa.String(length=1023), nullable=True),
    )
    op.execute("UPDATE ml_management_model SET cell_metadata_uri = 'gs://unknown'")
    op.alter_column("ml_management_model", "cell_metadata_uri", existing_type=sa.String(length=1023), nullable=False)


def downgrade() -> None:
    op.drop_column("ml_management_model", "cell_metadata_uri")
