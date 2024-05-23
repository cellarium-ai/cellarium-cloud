import os

MACHINE_SPEC_FILE_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../..", "settings/.pipelines-metafiles/machine_specs_tmp_data.json")
)
PCA_TRAIN_COMPONENT_NAME: str = "pca_train"
PCA_EMBED_COMPONENT_NAME: str = "pca_embed"
PCA_REGISTRY_COMPONENT_NAME: str = "pca_registry"
PCA_INDEX_CREATE_COMPONENT_NAME: str = "pca_index_create"
PCA_RESIZE_AND_SAVE_COMPONENT_NAME: str = "pca_resize_and_save"
MEAN_VAR_STD_COMPONENT_NAME: str = "mean_var_std_train"
TDIGEST_COMPONENT_NAME: str = "tdigest_train"
TDIGEST_FILTER_FEATURES_COMPONENT_NAME: str = "tdigest_filter_features"
LOGISTIC_REGRESSION_TRAIN_COMPONENT_NAME: str = "logistic_regression_train"
BENCHMARKING_COMPONENT_NAME: str = "benchmarking"
BQ_OPS_CREATE_AVRO_FILES_COMPONENT_NAME: str = "bq_ops_create_avro_files"
BQ_OPS_INGEST_DATA_COMPONENT_NAME: str = "bq_ops_ingest_data"
BQ_OPS_PRECALCULATE_FIELDS_COMPONENT_NAME: str = "bq_ops_precalculate_fields"
BQ_OPS_PREPARE_EXTRACT_COMPONENT_NAME: str = "bq_ops_prepare_extract"
BQ_OPS_EXTRACT_COMPONENT_NAME: str = "bq_ops_extract"

NETWORK_NAME: str = "projects/350868384795/global/networks/ai-matching"
DEFAULT_PIPELINE_LOCATION: str = "us-central1"
