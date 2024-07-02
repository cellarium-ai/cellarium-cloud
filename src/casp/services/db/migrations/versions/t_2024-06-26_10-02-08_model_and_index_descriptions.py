"""model_and_index_descriptions

Revision ID: 05b03bfc2e69
Revises: 0045d3785268
Create Date: 2024-06-26 10:02:08.149516

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "05b03bfc2e69"
down_revision = "0045d3785268"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("ml_management_model", sa.Column("description", sa.String(511), nullable=True))
    op.add_column("ml_management_matchingengineindex", sa.Column("description", sa.String(511), nullable=True))


def downgrade() -> None:
    op.drop_column("ml_management_model", "description")
    op.drop_column("ml_management_matchingengineindex", "description")
