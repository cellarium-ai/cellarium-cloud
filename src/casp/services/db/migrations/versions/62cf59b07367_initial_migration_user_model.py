"""Initial Migration: User model

Revision ID: 62cf59b07367
Revises:
Create Date: 2023-02-06 12:21:00.028094

"""

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "62cf59b07367"
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table(
        "user",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("email", sa.String(length=255), nullable=False),
        sa.Column("username", sa.String(length=255), nullable=False),
        sa.Column("active", sa.Boolean(), nullable=False),
        sa.Column("requests_processed", sa.Integer(), nullable=False),
        sa.Column("cells_processed", sa.Integer(), nullable=False),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("email"),
    )
    # ### end Alembic commands ###


# Uncomment this code if you need to downgrade the migration
# def downgrade() -> None:
#     # ### commands auto generated by Alembic - please adjust! ###
#     op.drop_table('user')
#     # ### end Alembic commands ###
