"""user_activity

Revision ID: 310657f0e46c
Revises: 0009
Create Date: 2024-03-27 17:00:17.459888

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "310657f0e46c"
down_revision = "a2e0ae902a06"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "user_activity",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("user_id", sa.Integer(), nullable=True),
        sa.Column("cell_count", sa.Integer(), nullable=False),
        sa.Column("model_name", sa.String(length=255), nullable=False),
        sa.Column("method", sa.String(length=255), nullable=True),
        sa.Column("finished_time", sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(
            ["user_id"],
            ["user.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        sa.Index("ix_user_activity_user_id", "user_id", unique=False),
    )


def downgrade() -> None:
    op.drop_table("user_activity")
