from casp.services import settings

KNN_QUERY_TEMPLATE_DIR = f"{settings.BQ_SQL_TEMPLATES_DIR}/knn_query_match"
MATCH_METADATA_SQL_TMPL_PATH = f"{KNN_QUERY_TEMPLATE_DIR}/match_metadata.sql.mako"
MATCH_METADATA_DEV_DETAILS_SQL_TMPL_PATH = f"{KNN_QUERY_TEMPLATE_DIR}/match_metadata_dev_details.sql.mako"
