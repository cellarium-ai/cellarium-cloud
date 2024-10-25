import json
import os
import typing as t

import dotenv
from pydantic_settings import BaseSettings

# Directory that contains the services package
SERVICES_DIR = os.path.dirname(os.path.abspath(__file__))
# Directory that contains the CAS package's contents
CAS_DIR = os.path.dirname(SERVICES_DIR)
# Directory that contains the CAS package and its settings
ROOT_DIR = os.path.dirname(CAS_DIR)


dotenv.load_dotenv(dotenv_path=f"{ROOT_DIR}/settings/.env")

ENV_TYPE = os.environ.get("ENVIRONMENT")


class AllEnvSettings(BaseSettings):
    # General
    GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS", "{}"))
    ENVIRONMENT: str = ENV_TYPE
    APP_VERSION: str = "1.7.0"
    APP_ROOT: str = ROOT_DIR
    DEFAULT_FEATURE_SCHEMA: str = "refdata-gex-GRCh38-2020-A"
    PROJECT_BUCKET_NAME: t.Optional[str] = os.environ.get("PROJECT_BUCKET_NAME")
    SERVICES_DIR: str = SERVICES_DIR
    DEFAULT_SERVICE_HOST: str = "0.0.0.0"
    DEFAULT_SERVICE_PORT: int = 8000
    API_SERVICE_PORT: int = DEFAULT_SERVICE_PORT
    LOG_LEVEL: str = "info"
    LOG_CONFIG: str = f"{SERVICES_DIR}/log_config.yaml"
    LOG_AS_JSON: bool = os.environ.get("LOG_AS_JSON", True)

    # Sentry
    SENTRY_DSN: t.Optional[str] = os.environ.get("SENTRY_DSN")
    SENTRY_ENABLE_TRACING: bool = True
    SENTRY_PROFILES_SAMPLE_RATE: float = 1.0
    SENTRY_TRACES_SAMPLE_RATE: float = 1.0
    SENTRY_ENVIRONMENT: t.Optional[str] = os.environ.get("SENTRY_ENVIRONMENT")
    # Model Training
    NEPTUNE_API_KEY: t.Optional[str] = os.environ.get("NEPTUNE_API_KEY")
    # API
    AIOHTTP_CLIENT_TOTAL_TIMEOUT_SECONDS: int = 650  # 350 seconds
    AIOHTTP_CLIENT_READ_TIMEOUT_SECONDS: int = 600  # 300 seconds
    DEFAULT_SCHEMA_NAME: str = "refdata-gex-GRCh38-2020-A"
    DEFAULT_MODEL_BQ_DATASET_NAME: str = "cas_50m_dataset"
    API_REQUEST_TEMP_TABLE_DATASET: str = "dsp-cell-annotation-service.cellarium_api_temp_tables"
    API_REQUEST_TEMP_TABLE_DATASET_EXPIRATION: int = 10  # 10 minutes
    KNN_SEARCH_NUM_MATCHES_DEFAULT: int = 100
    ITEMS_PER_USER: int = 50
    GET_MATCHES_CHUNK_SIZE: int = 5
    GET_MATCHES_MAX_RETRIES: int = 10
    GET_MATCHES_RETRY_BACKOFF_MULTIPLIER: int = 2
    GET_MATCHES_RETRY_BACKOFF_MIN: int = 0
    GET_MATCHES_RETRY_BACKOFF_MAX: int = 30
    MAX_CELL_IDS_PER_QUERY: int = 20_000  # Maximum number of cell IDs that can be queried at once
    # Consensus Engine
    GCS_CELL_ONTOLOGY_RESOURCE_FILE: str = (
        f"gs://{PROJECT_BUCKET_NAME}/consensus_engine/cell_ontology_resources_schema_5.json"
    )
    # Auth
    JWT_HASHING_ALGORITHM: str = "HS256"
    JWT_DEFAULT_TOKEN_TTL: int = 60 * 60 * 24 * 180  # 180 days
    # Database
    DB_CONNECTION_POOL_SIZE: int = 50  # 50 connections
    DB_CONNECTION_POOL_MAX_OVERFLOW: int = 10  # 10 connections
    DB_CONNECTION_POOL_TIMEOUT: int = 40  # 40 seconds
    DB_CONNECTION_POOL_RECYCLE: int = 1800  # 30 minutes
    DB_CONNECTION_NAME: t.Optional[str] = os.environ.get("DB_CONNECTION_NAME")
    DB_NAME: t.Optional[str] = os.environ.get("DB_NAME")
    DB_PASSWORD: t.Optional[str] = os.environ.get("DB_PASSWORD")
    DB_USER: t.Optional[str] = os.environ.get("DB_USER")
    DB_INSTANCE_UNIX_SOCKET: t.Optional[str] = os.environ.get("DB_INSTANCE_UNIX_SOCKET")
    DB_PORT: t.Optional[str] = os.environ.get("DB_PORT")
    # If this is value is specified, it will be used instead of the unix socket
    DB_PRIVATE_IP: t.Optional[str] = os.environ.get("DB_PRIVATE_IP")
    # BigQuery
    BQ_SQL_TEMPLATES_DIR: str = f"{CAS_DIR}/datastore_manager/sql/templates"
    # Stage db connector through unix socket or private IP
    if DB_PRIVATE_IP is None:
        SQLALCHEMY_DATABASE_URI: str = (
            f"postgresql+pg8000://{DB_USER}:{DB_PASSWORD}@/{DB_NAME}?unix_sock={DB_INSTANCE_UNIX_SOCKET}/.s.PGSQL.5432"
        )
    else:
        SQLALCHEMY_DATABASE_URI: str = (
            f"postgresql+pg8000://{DB_USER}:{DB_PASSWORD}@{DB_PRIVATE_IP}:{DB_PORT}/{DB_NAME}"
        )
    # Admin
    SECRET_KEY: str = os.environ.get("FLASK_SECRET_KEY")
    SECURITY_PASSWORD_SALT: t.Optional[str] = os.environ.get("FLASK_SECURITY_PASSWORD_SALT")
    FLASK_ADMIN_SWATCH: str = "flatly"
    _FLASK_BASIC_AUTH_USERNAME: str = os.environ.get("FLASK_BASIC_AUTH_USERNAME", "")
    _FLASK_BASIC_AUTH_PASSWORD: str = os.environ.get("FLASK_BASIC_AUTH_PASSWORD", "")
    ADMIN_BASIC_AUTH_USER: t.Dict[str, str] = {
        "username": _FLASK_BASIC_AUTH_USERNAME,
        "password": _FLASK_BASIC_AUTH_PASSWORD,
    }
    DEBUG: bool = False
    # Email settings
    SENDGRID_API_KEY: str = os.environ.get("SENDGRID_API_KEY", "")
    FROM_ADDRESS: str = os.environ.get("FROM_ADDRESS", "Cellarium CAS <no-reply@cellarium.ai>")
    FEEDBACK_FORM_BASE_URL: t.Optional[str] = os.environ.get("FEEDBACK_FORM_BASE_URL")
    DB_LOG_QUERIES: bool = os.environ.get("DB_LOG_QUERIES", False)

    # Quota settings
    DEFAULT_LIFETIME_QUOTA: int = 50000
    INCREASED_LIFETIME_QUOTA: int = 200000


class DevSettings(AllEnvSettings):
    # General
    debug: bool = True


class ProductionSettings(AllEnvSettings):
    pass


class LocalSettings(AllEnvSettings):
    # General
    debug: bool = True
    API_SERVICE_PORT: int = 8000
    # Database
    DB_HOST: str = os.environ.get("DB_HOST")
    DB_PORT: str = os.environ.get("DB_PORT")
    DB_NAME: str = os.environ.get("DB_NAME")
    DB_PASSWORD: str = os.environ.get("DB_PASSWORD")
    DB_USER: str = os.environ.get("DB_USER")
    SQLALCHEMY_DATABASE_URI: str = f"postgresql+pg8000://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
    LOG_AS_JSON: bool = os.environ.get("LOG_AS_JSON", False)


class TestSettings(AllEnvSettings):
    SECRET_KEY: str = "test"
