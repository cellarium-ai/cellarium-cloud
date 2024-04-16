import logging
import typing as t

from google.auth import default
from google.cloud import storage
from google.cloud.storage import transfer_manager
from google.oauth2.service_account import Credentials

from casp.services import settings

logger = logging.getLogger(__name__)


def get_google_service_credentials() -> t.Tuple[Credentials, str]:
    if settings.GOOGLE_ACCOUNT_CREDENTIALS != {}:
        logger.warning(
            "Google auth performed using the GOOGLE_ACCOUNT_CREDENTIALS env setting. "
            + "GOOGLE_APPLICATION_CREDENTIALS env variable or native cloud identity should be used instead."
        )
        credentials = Credentials.from_service_account_info(
            info=settings.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
        )
        return credentials, settings.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")

    credentials, project_id = default()
    if "get_project_id" in project_id and project_id != credentials.get_project_id():
        logger.warning(
            (
                f"Active project {project_id} is different than the default credential's project {credentials.project_id}.",
                f"Using {credentials.project_id} to auth to Google service",
            )
        )
        project_id = credentials.project_id
    logger.info(f"Using Google Cloud project {project_id}")
    return credentials, project_id


def download_file_from_bucket(bucket_name: str, source_blob_name: str, destination_file_name: str) -> None:
    """
    Download a file from GCS

    :param bucket_name: Bucket name in Google Cloud Storage
    :param source_blob_name: Name of the source blob in Google Cloud Storage Bucket
    :param destination_file_name: Local file name where to save the blob
    """
    storage_client = storage.Client()
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
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    blob = bucket.blob(blob_name=source_blob_name)
    blob.download_to_file(file)


def upload_file_to_bucket(local_file_name: str, bucket: str, blob_name: str) -> None:
    """
    Upload local file to a remote GCP bucket

    :param local_file_name: Local file name to upload
    :param bucket: Bucket name in Google Cloud Storage
    :param blob_name: Name of the destination blob in Bucket
    """
    client = storage.Client()
    bucket = client.get_bucket(bucket)
    blob = bucket.blob(blob_name)
    blob.upload_from_filename(local_file_name)


def upload_file_to_bucket_from_memory(bucket_name: str, contents: str, destination_blob_name: str):
    """
    Upload a file to the bucket

    :param bucket_name: Bucket name in GCS
    :param contents: Content which has to be written to blob
    :param destination_blob_name: Final destination in the bucket
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_string(contents)


def upload_many_blobs_with_transfer_manager(bucket_name, file_paths, prefix, workers=8):
    """Upload every file in a list to a bucket, concurrently in a process pool.

    Each blob name is derived from the filename, not including the
    `source_directory` parameter. For complete control of the blob name for each
    file (and other aspects of individual blob metadata), use
    transfer_manager.upload_many() instead.
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    filenames_blob_pairs = []

    for file_path in file_paths:
        file_name = file_path.split("/")[-1]
        blob = storage.Blob(bucket=bucket, name=f"{prefix}/{file_name}")
        filenames_blob_pairs.append((file_path, blob))

    results = transfer_manager.upload_many(file_blob_pairs=filenames_blob_pairs)

    number_of_uploaded = 0
    for name, result in zip(file_paths, results):
        # The results list is either `None` or an exception for each filename in
        # the input list, in order.

        if isinstance(result, Exception):
            print("Failed to upload {} due to exception: {}".format(name, result))
        else:
            number_of_uploaded += 1

    print(f"Successfully uploaded {number_of_uploaded} files")


def list_blobs(bucket_name: str, prefix: t.Optional[str] = None) -> t.List[storage.Blob]:
    """
    List all the blobs in the bucket

    :param bucket_name: Bucket name in Google Cloud Storage
    :param prefix: list all blobs with a prefix (a.k.a. "directory")

    :return: List of blobs with the provided prefix
    """
    storage_client = storage.Client()
    return storage_client.list_blobs(bucket_name, prefix=prefix)


def delete_folder_from_bucket(bucket_name: str, folder_name: str) -> None:
    """
    Delete a list of blobs with the same prefix

    :param bucket_name: Bucket name in Google Cloud Storage
    :param folder_name: Prefix which is used to mimic a folder. All files with the same prefixes will be removed
    """
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    blobs = list(bucket.list_blobs(prefix=folder_name))
    bucket.delete_blobs(blobs)
    print(f"Folder {folder_name} deleted")
