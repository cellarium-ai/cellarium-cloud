"""Add indexes to cell model for faster queries

Revision ID: 805c8eaafdf8
Revises: 86eda43c427a
Create Date: 2024-04-05 08:58:57.018405

"""

from alembic import op

# revision identifiers, used by Alembic.
revision = "805c8eaafdf8"
down_revision = "86eda43c427a"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_index(op.f("ix_cells_cellinfo_assay"), "cells_cellinfo", ["assay"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_cas_ingest_id"), "cells_cellinfo", ["cas_ingest_id"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_cell_type"), "cells_cellinfo", ["cell_type"], unique=False)
    op.create_index(
        op.f("ix_cells_cellinfo_cell_type_ontology_term_id"),
        "cells_cellinfo",
        ["cell_type_ontology_term_id"],
        unique=False,
    )
    op.create_index(op.f("ix_cells_cellinfo_development_stage"), "cells_cellinfo", ["development_stage"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_disease"), "cells_cellinfo", ["disease"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_is_primary_data"), "cells_cellinfo", ["is_primary_data"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_organism"), "cells_cellinfo", ["organism"], unique=False)
    op.create_index(
        op.f("ix_cells_cellinfo_self_reported_ethnicity"), "cells_cellinfo", ["self_reported_ethnicity"], unique=False
    )
    op.create_index(op.f("ix_cells_cellinfo_sex"), "cells_cellinfo", ["sex"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_suspension_type"), "cells_cellinfo", ["suspension_type"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_tissue"), "cells_cellinfo", ["tissue"], unique=False)
    op.create_index(op.f("ix_cells_cellinfo_total_mrna_umis"), "cells_cellinfo", ["total_mrna_umis"], unique=False)


def downgrade() -> None:
    op.drop_index(op.f("ix_cells_cellinfo_total_mrna_umis"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_tissue"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_suspension_type"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_sex"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_self_reported_ethnicity"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_organism"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_is_primary_data"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_disease"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_development_stage"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_cell_type_ontology_term_id"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_cell_type"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_cas_ingest_id"), table_name="cells_cellinfo")
    op.drop_index(op.f("ix_cells_cellinfo_assay"), table_name="cells_cellinfo")
