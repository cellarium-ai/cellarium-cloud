import json
import os
import typing as t

import dotenv
from pydantic import BaseSettings

dotenv.load_dotenv(dotenv_path="casp/admin/.env")


class Settings(BaseSettings):
    GOOGLE_ACCOUNT_CREDENTIALS: t.Dict = json.loads(os.environ.get("GOOGLE_SERVICE_ACCOUNT_CREDENTIALS"))
    SECRET_KEY: str = os.environ.get("FLASK_SECRET_KEY")
    SECURITY_PASSWORD_SALT = os.environ.get("FLASK_SECURITY_PASSWORD_SALT")
    _DB_HOST: str = os.environ.get("DB_HOST")
    _DB_PORT: str = os.environ.get("DB_PORT")
    _DB_NAME: str = os.environ.get("DB_NAME")
    _DB_PASSWORD: str = os.environ.get("DB_PASSWORD")
    _DB_USER: str = os.environ.get("DB_USER")
    SQLALCHEMY_DATABASE_URI: str = f"postgresql://{_DB_USER}:{_DB_PASSWORD}@{_DB_HOST}:{_DB_PORT}/{_DB_NAME}"
    FLASK_ADMIN_SWATCH: str = "flatly"
    FLASK_BASIC_AUTH_USERNAME: str = os.environ.get("FLASK_BASIC_AUTH_USERNAME")
    FLASK_BASIC_AUTH_PASSWORD: str = os.environ.get("FLASK_BASIC_AUTH_PASSWORD")
    DEBUG: bool = False if os.environ.get("ENVIRONMENT_TYPE") == "production" else True


settings = Settings()
