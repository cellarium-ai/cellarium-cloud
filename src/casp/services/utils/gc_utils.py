import typing as t

from google.cloud import storage
from google.oauth2.service_account import Credentials

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


def write_to_file_from_bucket(bucket_name: str, source_blob_name: str, file: t.BinaryIO) -> None:
    """
    Download the content of a file in GCS and write it to a file-like object

    :param bucket_name: Bucket name in Google Cloud Storage
    :param source_blob_name: Name of the source blob in Google Cloud Storage Bucket
    :param file: File to write from the blob
    """
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name)

    blob = bucket.blob(blob_name=source_blob_name)
    blob.download_to_file(file)


def upload_file_to_bucket(local_file_name: str, bucket: str, blob_name: str) -> None:
    credentials, project_id = get_google_service_credentials()
    client = storage.Client(credentials=credentials, project=project_id)
    bucket = client.get_bucket(bucket)
    blob = bucket.blob(blob_name)
    blob.upload_from_filename(local_file_name)


def list_blobs(bucket_name: str, prefix: t.Optional[str] = None) -> storage.Blob:
    """Lists all the blobs in the bucket"""
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(credentials=credentials, project=project_id)
    return storage_client.list_blobs(bucket_name, prefix=prefix)


def delete_folder_from_bucket(bucket_name: str, folder_name: str) -> None:
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(credentials=credentials, project=project_id)
    bucket = storage_client.get_bucket(bucket_name)
    blobs = list(bucket.list_blobs(prefix=folder_name))
    bucket.delete_blobs(blobs)
    print(f"Folder {folder_name} deleted")
