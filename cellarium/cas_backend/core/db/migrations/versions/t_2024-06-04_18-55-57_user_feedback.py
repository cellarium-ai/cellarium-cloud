"""user feedback

Revision ID: dbee57c2946d
Revises: 48b338e2e8e0
Create Date: 2024-06-04 18:55:57.551963

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "dbee57c2946d"
down_revision = "48b338e2e8e0"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("users_user", sa.Column("ask_for_feedback", sa.Boolean(), nullable=True))
    # Backfill the new column's data with default value (e.g. is_grpc = True)
    op.execute("update users_user set ask_for_feedback = true")
    # Change the new column to not nullable
    op.alter_column("users_user", "ask_for_feedback", existing_type=sa.BOOLEAN(), nullable=False)


def downgrade() -> None:
    op.drop_column("users_user", "ask_for_feedback")
