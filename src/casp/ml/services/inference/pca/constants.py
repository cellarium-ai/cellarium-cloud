import json
import os
import typing as t

import dotenv

dotenv.load_dotenv(dotenv_path="casp/ml/inference/pca/.env")
BUCKET_NAME: str = "fedor-test-bucket"
BLOB_NAME: str = "models/dump_manager_002.pickle"
GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
