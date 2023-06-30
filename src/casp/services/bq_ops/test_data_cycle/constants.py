GCS_INPUT_PATHS_HOMO_SAPIENS = [
    "test_input/tg66rtgh7y-test-source-file-1.h5ad",
    "test_input/yt55tgy4o-test-source-file-3.h5ad",
    "test_input/sx9871has78os-test-source-file-4.h5ad",
    "test_input/polloj3d39-test-source-file-5.h5ad",
]
GCS_INPUT_PATHS_MUS_MUS = ["test_input/tes2dc49-test-source-file-mus-mus-6.h5ad"]
GCS_INPUT_FILE_PATHS_ALL = [
    *GCS_INPUT_PATHS_HOMO_SAPIENS,
    *GCS_INPUT_PATHS_MUS_MUS,
]
DATASET_NAME = "cas_test_dataset"
HOMO_SAPIENS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A"
MUS_MUS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-mm10-2020-A"
ALL_CELLS_EXTRACT_TABLE_PREFIX = "test_extract_all"
HOMO_SAPIENS_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap"
HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap_5k"
FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX = "test_extract_filter_by_datasets"
MUS_MUS_EXTRACT_TABLE_PREFIX = "test_extract_mus_mus"
GCS_BUCKET_NAME = "dsp-cell-annotation-service"
GCS_STAGE_DIR = "test_tmp_dir"
