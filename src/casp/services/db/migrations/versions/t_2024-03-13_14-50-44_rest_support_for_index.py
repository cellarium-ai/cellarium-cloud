"""rest_support_for_index

Revision ID: 2633d266e26d
Revises: 0009
Create Date: 2024-03-13 14:50:44.326947

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "2633d266e26d"
down_revision = "0009"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("cas_matching_engine_index", sa.Column("is_grpc", sa.Boolean(), nullable=True))
    # Backfill the new column's data with default value (e.g. is_grpc = True)
    op.execute("update cas_matching_engine_index set is_grpc = true")
    # Change the new column to not nullable
    op.alter_column("cas_matching_engine_index", "is_grpc", existing_type=sa.BOOLEAN(), nullable=False)

    op.add_column("cas_matching_engine_index", sa.Column("api_endpoint", sa.String(length=255), nullable=True))

    # Retroactively apply.  These should be non-nullalble but must have been missed in the original migration.
    op.alter_column("cas_matching_engine_index", "index_name", existing_type=sa.VARCHAR(length=255), nullable=False)
    op.alter_column("cas_model", "model_name", existing_type=sa.VARCHAR(length=255), nullable=False)


def downgrade() -> None:
    op.alter_column("cas_model", "model_name", existing_type=sa.VARCHAR(length=255), nullable=True)
    op.alter_column("cas_matching_engine_index", "index_name", existing_type=sa.VARCHAR(length=255), nullable=True)
    op.drop_column("cas_matching_engine_index", "api_endpoint")
    op.drop_column("cas_matching_engine_index", "is_grpc")
