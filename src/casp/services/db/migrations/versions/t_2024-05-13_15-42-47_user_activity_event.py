"""user_activity_event

Revision ID: 9f6eb1c6a5a0
Revises: 22dce23f1475
Create Date: 2024-05-13 15:42:47.232019

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "9f6eb1c6a5a0"
down_revision = "22dce23f1475"
branch_labels = None
depends_on = None


def upgrade() -> None:
    user_activity_event = sa.Enum("SUCCEEDED", "FAILED", "STARTED", name="useractivityevent")
    user_activity_event.create(op.get_bind())

    op.add_column("users_useractivity", sa.Column("event", user_activity_event))

    op.execute("update users_useractivity set event = 'SUCCEEDED'")
    op.alter_column("users_useractivity", "event", existing_type=user_activity_event, nullable=False)

    op.create_index(op.f("ix_users_useractivity_event"), "users_useractivity", ["event"], unique=False)

    op.add_column("users_useractivity", sa.Column("request_id", sa.String(255), nullable=True))


def downgrade() -> None:
    op.drop_column("users_useractivity", "request_id")

    op.drop_index(op.f("ix_users_useractivity_event"), table_name="users_useractivity")

    op.drop_column("users_useractivity", "event")

    user_activity_event = sa.Enum("SUCCEEDED", "FAILED", "STARTED", name="useractivityevent")
    user_activity_event.drop(op.get_bind())
