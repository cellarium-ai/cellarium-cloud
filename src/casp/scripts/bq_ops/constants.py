from casp.services import settings

MRNA_FEATURE_BIOTYPE_NAME = "gene"
CAS_CELL_INFO_REQUIRED_COLUMNS = ["c.cas_cell_index", "c.cas_ingest_id"]
# Directories
SQL_TEMPLATES_DIR = f"{settings.SCRIPTS_DIR}/bq_ops/sql_templates"
PREPARE_CURRICULUM_SQL_DIR = f"{SQL_TEMPLATES_DIR}/prepare_curriculum"
EXTRACT_CURRICULUM_SQL_DIR = f"{SQL_TEMPLATES_DIR}/extract_curriculum"
# SQL Templates
PREPARE_CELL_INFO_RAND_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_cell_info_randomized.sql.mako"
PREPARE_CELL_INFO_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_cell_info.sql.mako"
DROP_PREPARE_CI_RAND_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/drop_prepare_cell_info_randomized.sql.mako"
PREPARE_FEATURE_SUMMARY_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_feature_summary.sql.mako"
PREPARE_CATEGORICAL_VARIABLE_SQL_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_categorical_variable.sql.mako"
GET_CELLS_IN_BIN_RANGE_SQL_DIR = f"{EXTRACT_CURRICULUM_SQL_DIR}/get_cells_in_bin_range.sql.mako"
# Metadata File name formats and strings
SHARED_META_DIR_NAME = "shared_meta"
MEASURED_GENES_INFO_FILE_NAME = "measured_genes_info.csv"
CATEGORICAL_COLUMNS_META_FILE_NAME = "categorical_columns_meta.pickle"
