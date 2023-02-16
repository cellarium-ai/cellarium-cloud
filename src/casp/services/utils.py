import typing as t

from google.oauth2.service_account import Credentials
from google.cloud import storage

from casp.services import settings


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")


def download_file_from_bucket(bucket_name: str, source_blob_name: str, destination_file_name: str) -> None:
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name)

    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)
