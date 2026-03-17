"""JSON -> JSONB

Revision ID: 429a028f2c30
Revises: 05cf37c7eda9
Create Date: 2024-04-11 12:49:54.739691

"""

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = "429a028f2c30"
down_revision = "05cf37c7eda9"
branch_labels = None
depends_on = None


def upgrade() -> None:
    with op.batch_alter_table("cells_cellinfo", schema=None) as batch_op:
        batch_op.alter_column(
            "obs_metadata_extra", type_=JSONB, existing_type=sa.JSON, postgresql_using="obs_metadata_extra::jsonb"
        )

    with op.batch_alter_table("cells_ingestinfo", schema=None) as batch_op:
        batch_op.alter_column(
            "uns_metadata", type_=JSONB, existing_type=sa.JSON, postgresql_using="uns_metadata::jsonb"
        )

    with op.batch_alter_table("cells_featureinfo", schema=None) as batch_op:
        batch_op.alter_column(
            "var_metadata_extra", type_=JSONB, existing_type=sa.JSON, postgresql_using="var_metadata_extra::jsonb"
        )


def downgrade():
    with op.batch_alter_table("cells_cellinfo", schema=None) as batch_op:
        batch_op.alter_column(
            "obs_metadata_extra", type_=sa.JSON, existing_type=JSONB, postgresql_using="obs_metadata_extra::json"
        )
    with op.batch_alter_table("cells_ingestinfo", schema=None) as batch_op:
        batch_op.alter_column("uns_metadata", type_=sa.JSON, existing_type=JSONB, postgresql_using="uns_metadata::json")
    with op.batch_alter_table("cells_featureinfo", schema=None) as batch_op:
        batch_op.alter_column(
            "var_metadata_extra", type_=sa.JSON, existing_type=JSONB, postgresql_using="var_metadata_extra::json"
        )
