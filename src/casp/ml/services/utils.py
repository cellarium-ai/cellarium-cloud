import typing as t

from google.oauth2.service_account import Credentials

from casp.ml.services import settings


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")
