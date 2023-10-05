"""Cas Model BQ table and more

Revision ID: 75ad185b37c1
Revises: 3518a8f30bcd
Create Date: 2023-06-23 12:12:21.111349

"""
import sqlalchemy as sa
from alembic import op

from casp.services import settings

# revision identifiers, used by Alembic.
revision = "75ad185b37c1"
down_revision = "3518a8f30bcd"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_column("cas_matching_engine_index", "admin_use_only")
    op.add_column("cas_model", sa.Column("schema_name", sa.String(length=255), nullable=True))
    op.add_column("cas_model", sa.Column("bq_cell_info_table_fqn", sa.String(length=255), nullable=True))
    op.add_column("cas_model", sa.Column("bq_temp_table_dataset", sa.String(length=255), nullable=True))
    # Setting up default values:
    op.execute(
        f"""
        UPDATE cas_model
            SET     schema_name = '{settings.DEFAULT_SCHEMA_NAME}',
                    bq_cell_info_table_fqn = '{settings.DEFAULT_MODEL_CELL_INFO_TABLE_FQN}',
                    bq_temp_table_dataset = '{settings.DEFAULT_MODEL_BQ_TEMP_TABLE_DATASET}'
    """
    )
    op.alter_column("cas_model", "schema_name", nullable=False)
    op.alter_column("cas_model", "bq_cell_info_table_fqn", nullable=False)
    op.alter_column("cas_model", "bq_temp_table_dataset", nullable=False)
    # ### end Alembic commands ###


# def downgrade() -> None:
#     # ### commands auto generated by Alembic - please adjust! ###
#     op.drop_column("cas_model", "bq_temp_table_dataset")
#     op.drop_column("cas_model", "bq_cell_info_table_fqn")
#     op.drop_column("cas_model", "schema_name")
#     op.add_column(
#         "cas_matching_engine_index", sa.Column("admin_use_only", sa.BOOLEAN(), autoincrement=False, nullable=False)
#     )
#     # ### end Alembic commands ###


# %%