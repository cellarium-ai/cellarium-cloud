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


settings = Settings()
