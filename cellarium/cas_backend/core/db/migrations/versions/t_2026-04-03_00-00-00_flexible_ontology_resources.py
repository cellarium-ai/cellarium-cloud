"""flexible_ontology_resources

Revision ID: b2c3d4e5f6a7
Revises: a1b2c3d4e5f6
Create Date: 2026-04-03 00:00:00.000000

"""

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "b2c3d4e5f6a7"
down_revision = "a1b2c3d4e5f6"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # 1. Create ml_management_cell_info_metadata table
    op.create_table(
        "ml_management_cell_info_metadata",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("soma_dataframe_uri", sa.String(length=1023), nullable=False),
        sa.Column("created_date", sa.DateTime(), nullable=True),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("name"),
    )

    # 2. Create ml_management_ontological_column table
    op.create_table(
        "ml_management_ontological_column",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("cell_info_metadata_id", sa.Integer(), nullable=False),
        sa.Column("column_name", sa.String(length=255), nullable=False),
        sa.Column("ontology_resource_uri", sa.String(length=1023), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.ForeignKeyConstraint(
            ["cell_info_metadata_id"],
            ["ml_management_cell_info_metadata.id"],
        ),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("cell_info_metadata_id", "column_name"),
    )

    # 3. Add cell_info_metadata_id column (nullable initially for data migration)
    op.add_column(
        "ml_management_model",
        sa.Column("cell_info_metadata_id", sa.Integer(), nullable=True),
    )

    # 4. Migrate existing cell_metadata_uri data into CellInfoMetadata records
    op.execute("""
        INSERT INTO ml_management_cell_info_metadata (name, soma_dataframe_uri, created_date)
        SELECT model_name || ' metadata', cell_metadata_uri, NOW()
        FROM ml_management_model
    """)
    op.execute("""
        UPDATE ml_management_model m
        SET cell_info_metadata_id = cim.id
        FROM ml_management_cell_info_metadata cim
        WHERE cim.name = m.model_name || ' metadata'
    """)

    # 5. Add FK constraint and make column NOT NULL
    op.create_foreign_key(
        "fk_ml_management_model_cell_info_metadata",
        "ml_management_model",
        "ml_management_cell_info_metadata",
        ["cell_info_metadata_id"],
        ["id"],
    )
    op.alter_column("ml_management_model", "cell_info_metadata_id", existing_type=sa.Integer(), nullable=False)

    # 6. Drop the old cell_metadata_uri column
    op.drop_column("ml_management_model", "cell_metadata_uri")


def downgrade() -> None:
    # 1. Add back cell_metadata_uri (nullable initially for data migration)
    op.add_column(
        "ml_management_model",
        sa.Column("cell_metadata_uri", sa.String(length=1023), nullable=True),
    )

    # 2. Restore data from linked CellInfoMetadata
    op.execute("""
        UPDATE ml_management_model m
        SET cell_metadata_uri = cim.soma_dataframe_uri
        FROM ml_management_cell_info_metadata cim
        WHERE m.cell_info_metadata_id = cim.id
    """)

    # 3. Make NOT NULL
    op.alter_column("ml_management_model", "cell_metadata_uri", existing_type=sa.String(length=1023), nullable=False)

    # 4. Drop FK constraint and cell_info_metadata_id column
    op.drop_constraint("fk_ml_management_model_cell_info_metadata", "ml_management_model", type_="foreignkey")
    op.drop_column("ml_management_model", "cell_info_metadata_id")

    # 5. Drop child table before parent
    op.drop_table("ml_management_ontological_column")
    op.drop_table("ml_management_cell_info_metadata")
