"""increase default quota

Revision ID: 2e2ff7bcff8a
Revises: 51a4ad2bb22d
Create Date: 2024-07-22 17:14:21.590727

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "2e2ff7bcff8a"
down_revision = "51a4ad2bb22d"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.execute("update users_user set cell_quota = 100000 where cell_quota = 50000")
    op.alter_column("users_user", "cell_quota", existing_type=sa.INTEGER(), default=100000)


def downgrade() -> None:
    # Don't undo the quota increases
    op.alter_column("users_user", "cell_quota", existing_type=sa.INTEGER(), default=50000)
