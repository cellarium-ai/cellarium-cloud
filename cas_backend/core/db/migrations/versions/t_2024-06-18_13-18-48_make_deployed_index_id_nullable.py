"""Make deployed_index_id nullable

Revision ID: 0045d3785268
Revises: dbee57c2946d
Create Date: 2024-06-18 13:18:48.605205

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "0045d3785268"
down_revision = "dbee57c2946d"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.alter_column(
        "ml_management_matchingengineindex", "deployed_index_id", existing_type=sa.VARCHAR(length=255), nullable=True
    )


def downgrade() -> None:
    op.alter_column(
        "ml_management_matchingengineindex", "deployed_index_id", existing_type=sa.VARCHAR(length=255), nullable=False
    )
