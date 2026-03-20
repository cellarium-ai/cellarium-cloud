"""rename cell related models and tables

Revision ID: 05cf37c7eda9
Revises: 805c8eaafdf8
Create Date: 2024-04-11 11:01:21.648497

"""

from alembic import op

# revision identifiers, used by Alembic.
revision = "05cf37c7eda9"
down_revision = "805c8eaafdf8"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.rename_table("cells_cellfeature", "cells_featureinfo")
    op.rename_table("cells_cellingestinfo", "cells_ingestinfo")
    op.rename_table("cas_model", "ml_management_model")
    op.rename_table("cas_matching_engine_index", "ml_management_matchingengineindex")
    op.rename_table("user", "users_user")
    op.rename_table("user_activity", "users_useractivity")


def downgrade() -> None:
    op.rename_table("cells_featureinfo", "cells_cellfeature")
    op.rename_table("cells_ingestinfo", "cells_cellingestinfo")
    op.rename_table("ml_management_model", "cas_model")
    op.rename_table("ml_management_matchingengineindex", "cas_matching_engine_index")
    op.rename_table("users_user", "user")
    op.rename_table("users_useractivity", "user_activity")
