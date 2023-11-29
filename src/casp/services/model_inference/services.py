import pickle
import typing as t

import anndata
import torch
from google.cloud import storage

from casp.ml.dump_manager import DumpManager
from casp.services import settings, utils
from casp.services.model_inference import schemas
from casp.services.model_inference.data_managers import ModelInferenceDataManager


class ModelInferenceService:
    """Service for model inference."""

    def __init__(self):
        self.model_inference_dm = ModelInferenceDataManager()

    @staticmethod
    def get_dump_manager(dump_manager_location: str) -> DumpManager:
        """
        Get Model Dump Manager from Google Cloud Storage

        :param dump_manager_location: Location of model dump manager in Google Cloud Storage

        :return: DumpManager object
        """
        credentials, project_id = utils.get_google_service_credentials()
        storage_client = storage.Client(project=project_id, credentials=credentials)
        bucket = storage_client.bucket(bucket_name=settings.PROJECT_BUCKET_NAME)
        blob = bucket.blob(dump_manager_location)
        pickle_in = blob.download_as_string()
        return pickle.loads(pickle_in)

    @staticmethod
    def load_data(file: t.BinaryIO) -> t.Tuple[torch.Tensor, list]:
        """
        Load Sparse Matrix and ids from anndata file

        :param file: File of :class:`anndata.Anndata` object to load

        :return: Tuple of Sparse Matrix and ids
        """
        adata = anndata.read_h5ad(file)
        X = torch.Tensor(adata.X.todense().astype(int))
        obs_ids = adata.obs.index.values.astype(str).tolist()
        return X, obs_ids

    def embed_adata_file(self, file_to_embed: t.BinaryIO, model_name: str) -> schemas.ModelEmbeddings:
        """
        Embed adata file using a specific model using Cellarium-ML model and pytorch

        :param file_to_embed: File object of :class:`anndata.AnnData` object to embed.
        :param model_name: Model name to use for embedding.

        """
        # Get Model dump file
        model_info = self.model_inference_dm.get_model_by(model_name=model_name)
        dump_manager = self.get_dump_manager(model_info.model_file_path)
        model = dump_manager.model
        transform = dump_manager.transform
        # Get data and transform it
        X, obs_ids = self.load_data(file_to_embed)
        X = transform(X)
        # Embed data
        embeddings = model.transform(X).numpy()

        embeddings_b64 = utils.numpy_to_base64(embeddings)

        return schemas.ModelEmbeddings(obs_ids=obs_ids, embeddings_b64=embeddings_b64)
