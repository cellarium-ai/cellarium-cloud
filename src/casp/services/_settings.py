import json
import os
import typing as t

import dotenv
from pydantic import BaseSettings

dotenv.load_dotenv(dotenv_path="casp/services/.env")

ENV_TYPE = os.environ.get("ENVIRONMENT")


class AllEnvSettings(BaseSettings):
    # General
    GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS", "{}"))
    ENVIRONMENT: str = os.environ.get("ENVIRONMENT")
    # PCA Model serving
    PCA_MODEL_BUCKET_NAME: str = "fedor-test-bucket"
    PCA_MODEL_BLOB_NAME: str = "models/dump_manager_002.pickle"
    # Model Training
    NEPTUNE_API_KEY: str = os.environ.get("NEPTUNE_API_KEY")
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
    # Auth
    JWT_HASHING_ALGORITHM: str = "HS256"
    JWT_DEFAULT_TOKEN_TTL: int = 60 * 60 * 24 * 180  # 180 days
    # Database
    DB_NAME: str = os.environ.get("DB_NAME")
    DB_PASSWORD: str = os.environ.get("DB_PASSWORD")
    DB_USER: str = os.environ.get("DB_USER")
    DB_INSTANCE_UNIX_SOCKET: str = os.environ.get("DB_INSTANCE_UNIX_SOCKET")
    # Stage db connector through unix socket
    SQLALCHEMY_DATABASE_URI: str = (
        f"postgresql+pg8000://{DB_USER}:{DB_PASSWORD}@/{DB_NAME}?unix_sock={DB_INSTANCE_UNIX_SOCKET}/.s.PGSQL.5432"
    )
    # Admin
    SECRET_KEY: str = os.environ.get("FLASK_SECRET_KEY")
    SECURITY_PASSWORD_SALT: str = os.environ.get("FLASK_SECURITY_PASSWORD_SALT")
    FLASK_ADMIN_SWATCH: str = "flatly"
    _FLASK_BASIC_AUTH_USERNAME: str = os.environ.get("FLASK_BASIC_AUTH_USERNAME")
    _FLASK_BASIC_AUTH_PASSWORD: str = os.environ.get("FLASK_BASIC_AUTH_PASSWORD")
    ADMIN_BASIC_AUTH_USER: t.Dict[str, str] = {
        "username": _FLASK_BASIC_AUTH_USERNAME,
        "password": _FLASK_BASIC_AUTH_PASSWORD,
    }
    DEBUG: bool = False


class DevSettings(AllEnvSettings):
    # General
    debug = True


class ProductionSettings(AllEnvSettings):
    pass


class LocalSettings(AllEnvSettings):
    # General
    debug = True
    # Database
    DB_HOST: str = os.environ.get("DB_HOST")
    DB_PORT: str = os.environ.get("DB_PORT")
    DB_NAME: str = os.environ.get("DB_NAME")
    DB_PASSWORD: str = os.environ.get("DB_PASSWORD")
    DB_USER: str = os.environ.get("DB_USER")
    SQLALCHEMY_DATABASE_URI = f"postgresql+pg8000://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
