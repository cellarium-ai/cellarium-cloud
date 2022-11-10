import json
import os
import typing as t

import dotenv

dotenv.load_dotenv(dotenv_path="casp/ml/services/training/.env")
GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
