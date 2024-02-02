import typing as t

import anndata
import numpy
from cellarium.ml import CellariumModule
from cellarium.ml.utilities.data import AnnDataField, collate_fn
from jsonargparse import Namespace

from casp.services import settings, utils
from casp.services.db import models
from casp.services.model_inference import schemas
from casp.services.model_inference.data_managers import ModelInferenceDataManager


class ModelInferenceService:
    """Service for model inference."""

    def __init__(self):
        self.model_inference_dm = ModelInferenceDataManager()

    @staticmethod
    def _get_batch_dictionary_from_adata_file(
        adata_file: t.BinaryIO, cellarium_module: CellariumModule
    ) -> t.Dict[str, str]:
        """
        Get batch dictionary from adata file using cellarium module.

        :param adata_file: File object of :class:`anndata.AnnData` object to embed.
        :param cellarium_module: CellariumModule object

        :return: Batch dictionary
        """
        batch_keys = {}
        _batch_keys = cellarium_module._hparams["data"]["batch_keys"]
        _batch_keys["obs_names"] = Namespace(attr="obs_names")

        for key, value in vars(_batch_keys).items():
            batch_keys[key] = AnnDataField(**vars(value))

        adata = anndata.read_h5ad(adata_file)
        batch = {key: field(adata)[:] for key, field in batch_keys.items()}
        batch = collate_fn([batch])
        return batch

    @staticmethod
    def _get_model_checkpoint_path(model_checkpoint_file_path: str) -> str:
        """
        Get model checkpoint path from GCS. It expects a filepath that doesn't include bucket protocol prefix and
        bucket name.

        :param model_checkpoint_file_path: Model checkpoint file path in GCS

        :return: Model checkpoint file path with GCS protocol prefix
        """
        return f"gs://{settings.PROJECT_BUCKET_NAME}/{model_checkpoint_file_path}"

    def _get_model_checkpoint(self, model: models.CASModel) -> CellariumModule:
        """
        Get model checkpoint from either local or GCS and load it using CellariumModule.

        :param model: Cellarium Cloud model db object

        :return: CellariumModule object
        """
        model_checkpoint_path = self._get_model_checkpoint_path(model.model_file_path)
        return CellariumModule.load_from_checkpoint(model_checkpoint_path, map_location="cpu")

    def _get_output_from_model(
        self, model: models.CASModel, adata_file: t.BinaryIO
    ) -> t.Tuple[numpy.ndarray, t.List[str]]:
        """
        Get output from cellarium-ml model that predicts embeddings given an input adata file.

        :param model: Cellarium Cloud model db object
        :param adata_file: File object of :class:`anndata.AnnData` object to embed.

        :return: Tuple of embeddings and obs_ids.
        """
        cellarium_module = self._get_model_checkpoint(model)

        batch = self._get_batch_dictionary_from_adata_file(adata_file=adata_file, cellarium_module=cellarium_module)
        cellarium_output_dict = cellarium_module(batch)

        embeddings = cellarium_output_dict["x_ng"].numpy()
        obs_ids = cellarium_output_dict["obs_names"].tolist()
        return embeddings, obs_ids

    def embed_adata_file(self, file_to_embed: t.BinaryIO, model_name: str) -> schemas.ModelEmbeddings:
        """
        Embed adata file using a specific model using Cellarium-ML model and pytorch

        :param file_to_embed: File object of :class:`anndata.AnnData` object to embed.
        :param model_name: Model name to use for embedding.

        :return: ModelEmbeddings schema object.
        """
        model_info = self.model_inference_dm.get_model_by(model_name=model_name)

        embeddings, obs_ids = self._get_output_from_model(model=model_info, adata_file=file_to_embed)

        embeddings_b64 = utils.numpy_to_base64(embeddings)

        return schemas.ModelEmbeddings(obs_ids=obs_ids, embeddings_b64=embeddings_b64)
