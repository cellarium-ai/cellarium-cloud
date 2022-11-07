import io
import pickle
import typing as t

import anndata
import torch
from google.oauth2.service_account import Credentials
from google.cloud import storage

from casp.ml.dump_manager import DumpManager
from casp.ml.inference.pca import constants


def get_google_service_credentials() -> t.Tuple["Credentials", str]:
    credentials = Credentials.from_service_account_info(
        info=constants.GOOGLE_ACCOUNT_CREDENTIALS, scopes=None, default_scopes=None
    )
    return credentials, constants.GOOGLE_ACCOUNT_CREDENTIALS.get("project_id")


def load_data(file) -> t.Tuple[torch.Tensor, torch.Tensor]:
    adata = anndata.read_h5ad(io.BytesIO(file))
    X = torch.Tensor(adata.raw.X.todense().astype(int))
    db_ids = torch.Tensor(adata.obs_names.values.astype(int))
    return X, db_ids


def get_dump_manager() -> "DumpManager":
    credentials, project_id = get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=constants.BUCKET_NAME)
    blob = bucket.blob(constants.BLOB_NAME)
    # TODO: Make PickleMixin capable of downloading directly from buckets
    pickle_in = blob.download_as_string()
    # blob.download_to_filename("./dm.pickle")
    return pickle.loads(pickle_in)
    # return pickle.load(open("./dm.pickle", "rb"))
