from casp.services import settings

MRNA_FEATURE_BIOTYPE_NAME = "gene"
CAS_CELL_INFO_REQUIRED_COLUMNS = ["c.cas_cell_index", "c.cas_ingest_id"]
PREPARE_CURRICULUM_SQL_DIR = f"{settings.BQ_SQL_TEMPLATES_DIR}/prepare_curriculum"
EXTRACT_CURRICULUM_SQL_DIR = f"{settings.BQ_SQL_TEMPLATES_DIR}/extract_curriculum"
PREPARE_CELL_INFO_RAND_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_cell_info_randomized.sql.mako"
PREPARE_CELL_INFO_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_cell_info.sql.mako"
DROP_PREPARE_CI_RAND_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/drop_prepare_cell_info_randomized.sql.mako"
PREPARE_FEATURE_SUMMARY_TEMPLATE_DIR = f"{PREPARE_CURRICULUM_SQL_DIR}/prepare_feature_summary.sql.mako"
GET_CELLS_IN_BIN_RANGE_SQL_DIR = f"{EXTRACT_CURRICULUM_SQL_DIR}/get_cells_in_bin_range.sql.mako"
