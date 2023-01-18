import re

import anndata
from google.cloud import storage


def read_h5ad_gcs(
    filename: str,
) -> anndata.AnnData:
    """
    Read `.h5ad`-formatted hdf5 file from the Google Cloud Storage.

    Parameters
    ----------
    filename
        File path to the data file in Cloud Storage.
    """
    # parse bucket and blob names from the filename
    if filename.startswith("gs:"):
        filename = re.sub(r"^gs://?", "", filename)
        bucket_name, blob_name = filename.split("/", 1)

        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(blob_name)

        with blob.open("rb") as f:
            return anndata.read_h5ad(f)
    if filename.startswith("file:"):
        filename = re.sub(r"^file://?", "", filename)
    return anndata.read_h5ad(filename)
