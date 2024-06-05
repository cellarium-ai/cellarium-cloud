"""user keys

Revision ID: 48b338e2e8e0
Revises: 429a028f2c30
Create Date: 2024-05-16 09:49:58.443718

"""

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "48b338e2e8e0"
down_revision = "9f6eb1c6a5a0"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "users_userkey",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("key_locator", postgresql.UUID(), nullable=False),
        sa.Column("user_id", sa.Integer(), nullable=False),
        sa.Column("created_date", sa.DateTime(), nullable=False),
        sa.Column("active", sa.Boolean(), nullable=False),
        sa.Column("expires", sa.DateTime(), nullable=False),
        sa.Column("key_hash", sa.Text(), nullable=False),
        sa.ForeignKeyConstraint(
            ["user_id"],
            ["users_user.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("key_locator"),
    )


def downgrade() -> None:
    op.drop_table("users_userkey")
