"""add_memory_budget_to_vector_index

Revision ID: c3d4e5f6a7b8
Revises: b2c3d4e5f6a7
Create Date: 2026-04-15 00:00:00.000000

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "c3d4e5f6a7b8"
down_revision = "b2c3d4e5f6a7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("ml_management_vectorindex", sa.Column("memory_budget", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("ml_management_vectorindex", "memory_budget")
