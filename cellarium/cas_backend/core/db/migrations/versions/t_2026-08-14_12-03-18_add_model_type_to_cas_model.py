"""add_model_type_to_cas_model

Revision ID: 2b456929aacc
Revises: 6426fc752f63
Create Date: 2026-08-14 12:03:18.864384

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "2b456929aacc"
down_revision = "6426fc752f63"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("ml_management_model", sa.Column("model_type", sa.String(length=255), nullable=True))
    op.execute("update ml_management_model set model_type = 'representation'")
    op.alter_column("ml_management_model", "model_type", existing_type=sa.String(length=255), nullable=False)
    op.alter_column("ml_management_model", "embedding_dimension", existing_type=sa.INTEGER(), nullable=True)


def downgrade() -> None:
    op.execute("update ml_management_model set embedding_dimension = 0 where embedding_dimension is null")
    op.alter_column("ml_management_model", "embedding_dimension", existing_type=sa.INTEGER(), nullable=False)
    op.drop_column("ml_management_model", "model_type")
