GCS_INPUT_PATHS_HOMO_SAPIENS = [
    "dev/test_input/tg66rtgh7y__test-source-file-1.h5ad",
    "dev/test_input/yt55tgy4o__test-source-file-3.h5ad",
    "dev/test_input/sx9871has78os__test-source-file-4.h5ad",
    "dev/test_input/polloj3d39__test-source-file-5.h5ad",
]
GCS_INPUT_PATHS_MUS_MUS = ["dev/test_input/tes2dc49__test-source-file-mus-mus-6.h5ad"]
GCS_INPUT_FILE_PATHS_ALL = [
    *GCS_INPUT_PATHS_HOMO_SAPIENS,
    *GCS_INPUT_PATHS_MUS_MUS,
]
GCS_EXTRACT_FILE_NAME_PREFIX = "extract_{chunk_id}.h5ad"
DATASET_NAME = "cas_test_dataset"
HOMO_SAPIENS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A"
MUS_MUS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-mm10-2020-A"
ALL_CELLS_EXTRACT_TABLE_PREFIX = "test_extract_all"
HOMO_SAPIENS_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap"
MUS_MUS_EXTRACT_TABLE_PREFIX = "test_extract_mus_mus"
HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap_5k"
FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX = "test_extract_filter_by_datasets"
FILTER_BY_DISEASES_EXTRACT_TABLE_PREFIX = "test_extract_filter_by_diseases"
FILTERS_JSON_DIR = "tests/infrastructure/filters"
FILTER_HOMO_SAP_JSON_PATH = f"{FILTERS_JSON_DIR}/homo_sapiens.json"
FILTER_MUS_MUS_JSON_PATH = f"{FILTERS_JSON_DIR}/mus_musculus.json"
FILTER_HOMO_SAP_NO_CANCER_JSON_PATH = f"{FILTERS_JSON_DIR}/homo_sapiens_no_cancer.json"
FILTER_DATASET_FILENAME_JSON_PATH = f"{FILTERS_JSON_DIR}/filter_by_datasets.json"
ALL_CELLS_EXTRACT_BUCKET_PATH = f"curriculum/{ALL_CELLS_EXTRACT_TABLE_PREFIX}"
HOMO_SAPIENS_EXTRACT_BUCKET_PATH = f"curriculum/{HOMO_SAPIENS_EXTRACT_TABLE_PREFIX}"
MUS_MUS_EXTRACT_BUCKET_PATH = f"curriculum/{MUS_MUS_EXTRACT_TABLE_PREFIX}"
HOMO_SAPIENS_5k_EXTRACT_BUCKET_PATH = f"curriculum/{HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX}"
FILTER_BY_DATASET_EXTRACT_BUCKET_PATH = f"curriculum/{FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX}"
FILTER_BY_DISEASES_EXTRACT_BUCKET_PATH = f"curriculum/{FILTER_BY_DISEASES_EXTRACT_TABLE_PREFIX}"
GCS_BUCKET_NAME = "cellarium-file-system"
GCS_STAGE_DIR = "dev/test_tmp_dir"
PRECALCULATE_FIELDS = "total_mrna_umis"
OBS_COLUMNS_TO_INCLUDE = (
    "c.cell_type,c.total_mrna_umis,c.donor_id,c.assay,c.development_stage,"
    "c.disease,c.organism,c.sex,c.tissue,i.dataset_id"
)
OBS_COLUMNS_TO_INCLUDE_LIST = OBS_COLUMNS_TO_INCLUDE.split(",")
