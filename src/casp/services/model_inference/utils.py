import io
import pickle
import typing as t

import anndata
import torch
from google.cloud import storage

from casp.ml.dump_manager import DumpManager
from casp.services import settings, utils


def load_data(file) -> t.Tuple[torch.Tensor, torch.Tensor]:
    adata = anndata.read_h5ad(io.BytesIO(file))
    X = torch.Tensor(adata.X.todense().astype(int))
    db_ids = torch.Tensor(adata.obs_names.values.astype(int))
    return X, db_ids


def get_dump_manager(dump_manager_location: str) -> "DumpManager":
    credentials, project_id = utils.get_google_service_credentials()
    storage_client = storage.Client(project=project_id, credentials=credentials)
    bucket = storage_client.bucket(bucket_name=settings.PROJECT_BUCKET_NAME)
    blob = bucket.blob(dump_manager_location)
    pickle_in = blob.download_as_string()
    return pickle.loads(pickle_in)
