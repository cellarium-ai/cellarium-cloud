"""lifetime_quotas

Revision ID: c80bd73203d1
Revises: 2e2ff7bcff8a
Create Date: 2024-09-18 12:04:55.133331

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "c80bd73203d1"
down_revision = "2e2ff7bcff8a"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("users_user", sa.Column("lifetime_cell_quota", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("users_user", "lifetime_cell_quota")
