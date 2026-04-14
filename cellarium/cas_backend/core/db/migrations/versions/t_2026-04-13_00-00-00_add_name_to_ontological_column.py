"""add_name_to_ontological_column

Revision ID: c3d4e5f6a7b8
Revises: b2c3d4e5f6a7
Create Date: 2026-04-13 00:00:00.000000

"""

import uuid

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "c3d4e5f6a7b8"
down_revision = "b2c3d4e5f6a7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # 1. Add ontology_resource_name column as nullable first
    op.add_column(
        "ml_management_ontological_column",
        sa.Column("ontology_resource_name", sa.String(length=255), nullable=True),
    )

    # 2. Populate existing rows with a generated unique name
    connection = op.get_bind()
    rows = connection.execute(sa.text("SELECT id, column_name FROM ml_management_ontological_column")).fetchall()
    for row in rows:
        generated_name = f"{row[1]}_{uuid.uuid4().hex[:8]}"
        connection.execute(
            sa.text("UPDATE ml_management_ontological_column SET ontology_resource_name = :name WHERE id = :id"),
            {"name": generated_name, "id": row[0]},
        )

    # 3. Make column non-nullable and add unique constraint
    with op.batch_alter_table("ml_management_ontological_column") as batch_op:
        batch_op.alter_column("ontology_resource_name", existing_type=sa.String(length=255), nullable=False)
        batch_op.create_unique_constraint("uq_ontological_column_ontology_resource_name", ["ontology_resource_name"])


def downgrade() -> None:
    with op.batch_alter_table("ml_management_ontological_column") as batch_op:
        batch_op.drop_constraint("uq_ontological_column_ontology_resource_name", type_="unique")
        batch_op.drop_column("ontology_resource_name")
