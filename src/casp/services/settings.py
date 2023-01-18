import json
import os
import typing as t

import dotenv
from pydantic import BaseSettings

dotenv.load_dotenv(dotenv_path="casp/ml/services/.env")


class Settings(BaseSettings):
    # General
    GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
    # PCA Model serving:
    PCA_MODEL_BUCKET_NAME = "fedor-test-bucket"
    PCA_MODEL_BLOB_NAME = "models/dump_manager_002.pickle"
    # API
    SERVER_HOST: str = "0.0.0.0"
    SERVER_PORT: int = 8000
    MODEL_SERVER_URL: str = "https://casp-pca-serving-vi7nxpvk7a-uc.a.run.app/predict"
    KNN_SEARCH_ENDPOINT_ID: str = "projects/350868384795/locations/us-central1/indexEndpoints/2348891088464379904"
    KNN_SEARCH_DEPLOYED_INDEX_ID: str = "deployed_4m_casp_index_v1"
    KNN_SEARCH_NUM_MATCHES: int = 100
    BQ_CELL_INFO_TABLE_FQN: str = "dsp-cell-annotation-service.cas_4m_dataset.cas_cell_info"
    BQ_TEMP_TABLE_DATASET: str = "dsp-cell-annotation-service.cas_4m_dataset"
    ITEMS_PER_USER: int = 50


settings = Settings()
