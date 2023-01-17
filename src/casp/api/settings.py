import json
import os
import typing as t

import dotenv
from pydantic import BaseSettings

dotenv.load_dotenv(dotenv_path="casp/api/.env")


class Settings(BaseSettings):
    server_host: str = "0.0.0.0"
    server_port: int = 8000
    model_server_url: str = "https://casp-pca-serving-vi7nxpvk7a-uc.a.run.app/predict"
    knn_search_endpoint_id: str = "projects/350868384795/locations/us-central1/indexEndpoints/2348891088464379904"
    knn_search_deployed_index_id: str = "deployed_4m_casp_index_v1"
    knn_search_num_matches: int = 100
    bq_cell_info_table_fqn: str = "dsp-cell-annotation-service.cas_4m_dataset.cas_cell_info"
    bq_temp_table_dataset: str = "dsp-cell-annotation-service.cas_4m_dataset"
    items_per_user: int = 50
    GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
    AUTH_HEADER: str = "Bearer"


settings = Settings()
