import typing as t

from google.cloud import storage
from google.oauth2.service_account import Credentials

from casp.services import settings


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")


def upload_file_to_bucket(local_file_name: str, bucket: str, blob_name: str) -> None:
    credentials, project_id = get_google_service_credentials()
    client = storage.Client(credentials=credentials, project=project_id)
    bucket = client.get_bucket(bucket)
    blob = bucket.blob(blob_name)
    blob.upload_from_filename(local_file_name)
