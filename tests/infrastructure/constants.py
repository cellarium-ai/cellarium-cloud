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
GCS_EXTRACT_FILE_NAME_PREFIX = "extract_{chunk_id}.h5ad"
DATASET_NAME = "cas_test_dataset"
HOMO_SAPIENS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A"
MUS_MUS_GENE_SCHEMA = "dsp-cell-annotation-service.cas_reference_data.refdata-gex-mm10-2020-A"
ALL_CELLS_EXTRACT_TABLE_PREFIX = "test_extract_all"
HOMO_SAPIENS_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap"
HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX = "test_extract_homo_sap_5k"
FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX = "test_extract_filter_by_datasets"
FILTER_BY_DISEASES_EXTRACT_TABLE_PREFIX = "test_extract_filter_by_diseases"
NO_CANCER_DISEASES_FILTER = (
    "age related macular degeneration 7,influenza,basal laminar drusen,Barrett esophagus,chronic kidney disease,"
    "endocrine pancreas disorder,hypersensitivity pneumonitis,post-COVID-19 disorder,temporal lobe epilepsy,"
    "cystic fibrosis,arrhythmogenic right ventricular cardiomyopathy,non-specific interstitial pneumonia,"
    "hyperplastic polyp,gastritis,cataract,lymphangioleiomyomatosis,acute kidney failure,"
    "common variable immunodeficiency,systemic lupus erythematosus,benign prostatic hyperplasia,"
    "chronic rhinitis,dementia,dilated cardiomyopathy,interstitial lung disease,stomach disorder,"
    "respiratory system disorder,Down syndrome,pulmonary fibrosis,type 1 diabetes mellitus,myocardial infarction,"
    "Crohn disease,chronic obstructive pulmonary disease,pulmonary emphysema,trisomy 18,non-compaction cardiomyopathy,"
    "hydrosalpinx,tubular adenoma,lymphadenitis,normal,pneumonia,pulmonary sarcoidosis,Crohn ileitis,COVID-19,"
    "Alzheimer disease,pilocytic astrocytoma,type 2 diabetes mellitus"
)
DATASETS_FILTER = "tg66rtgh7y-test-source-file-1_ingest_info.avro,yt55tgy4o-test-source-file-3_ingest_info.avro"
MUS_MUS_EXTRACT_TABLE_PREFIX = "test_extract_mus_mus"
GCS_BUCKET_NAME = "dsp-cell-annotation-service"
GCS_STAGE_DIR = "test_tmp_dir"
PRECALCULATE_FIELDS = "total_mrna_umis"
OBS_COLUMNS_TO_INCLUDE = "cell_type,total_mrna_umis,donor_id,assay,development_stage,disease,organism,sex,tissue"
OBS_COLUMNS_TO_INCLUDE_LIST = OBS_COLUMNS_TO_INCLUDE.split(",")
