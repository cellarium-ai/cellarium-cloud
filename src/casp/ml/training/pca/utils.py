import typing as t

from google.oauth2.service_account import Credentials

from casp.ml.training.pca import constants


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=constants.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, constants.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")
