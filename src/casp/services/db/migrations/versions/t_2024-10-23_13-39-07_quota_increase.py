"""quota_increase

Revision ID: cdfe110d1068
Revises: c80bd73203d1
Create Date: 2024-10-23 13:39:07.609718

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "cdfe110d1068"
down_revision = "c80bd73203d1"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("users_user", sa.Column("quota_increased", sa.Boolean(), nullable=False, server_default="False"))
    op.execute("update users_user set quota_increased = true where lifetime_cell_quota = 200000")


def downgrade() -> None:
    op.drop_column("users_user", "quota_increased")
