"""add_vector_index_model

Revision ID: 532fd37ec55a
Revises: cdfe110d1068
Create Date: 2026-03-23 15:00:00.000000

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "532fd37ec55a"
down_revision = "cdfe110d1068"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "ml_management_vectorindex",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("index_name", sa.String(length=255), nullable=False),
        sa.Column("description", sa.String(length=511), nullable=True),
        sa.Column("embedding_dimension", sa.Integer(), nullable=False),
        sa.Column("num_neighbors", sa.Integer(), nullable=False),
        sa.Column("index_uri", sa.String(length=1023), nullable=False),
        sa.Column("index_type", sa.String(length=255), nullable=False),
        sa.Column("distance_metric", sa.String(length=255), nullable=False),
        sa.Column("nprobe", sa.Integer(), nullable=True),
        sa.Column("l_search", sa.Integer(), nullable=True),
        sa.Column("model_id", sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(["model_id"], ["ml_management_model.id"]),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("index_name"),
        sa.UniqueConstraint("model_id"),
    )


def downgrade() -> None:
    op.drop_table("ml_management_vectorindex")
