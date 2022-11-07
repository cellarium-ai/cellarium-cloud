import typing as t
import os
import json
import dotenv

dotenv.load_dotenv(dotenv_path="casp/ml/training/pca/.env")
BUCKET_NAME: str = "fedor-test-bucket"
GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
