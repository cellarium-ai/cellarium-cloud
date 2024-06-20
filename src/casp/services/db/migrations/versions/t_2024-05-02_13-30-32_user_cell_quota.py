"""user_cell_quota

Revision ID: 22dce23f1475
Revises: 429a028f2c30
Create Date: 2024-05-02 13:30:32.439909

"""

import datetime

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "22dce23f1475"
down_revision = "429a028f2c30"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("users_user", sa.Column("cell_quota", sa.Integer(), nullable=True, default=50000))
    op.execute("update users_user set cell_quota = 50000")
    op.alter_column("users_user", "cell_quota", existing_type=sa.INTEGER(), nullable=False)

    op.add_column(
        "users_user",
        sa.Column("created_at", sa.DateTime(), nullable=True, default=datetime.datetime.now(datetime.timezone.utc)),
    )
    op.execute("update users_user set created_at = now()")
    op.alter_column("users_user", "created_at", existing_type=sa.DateTime(), nullable=False)


def downgrade() -> None:
    op.drop_column("users_user", "cell_quota")
    op.drop_column("users_user", "created_at")
