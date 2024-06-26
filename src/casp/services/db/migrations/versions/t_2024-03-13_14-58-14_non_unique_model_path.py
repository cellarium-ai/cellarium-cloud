"""non unique model path

Revision ID: a2e0ae902a06
Revises: 2633d266e26d
Create Date: 2024-03-13 14:58:14.174652

"""

from alembic import op

# revision identifiers, used by Alembic.
revision = "a2e0ae902a06"
down_revision = "2633d266e26d"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.drop_constraint("cas_model_model_file_path_key", "cas_model", type_="unique")


def downgrade() -> None:
    op.create_unique_constraint("cas_model_model_file_path_key", "cas_model", ["model_file_path"])
