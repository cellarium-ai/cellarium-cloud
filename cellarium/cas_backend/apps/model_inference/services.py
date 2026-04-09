from functools import cache
from io import BytesIO
import typing as t

import anndata
import numpy as np
from smart_open import open

from cellarium.cas_backend.apps.model_inference import exceptions
from cellarium.cas_backend.core.config import settings
from cellarium.cas_backend.core.db import models
from cellarium.ml import CellariumAnnDataDataModule, CellariumModule
from cellarium.ml.utilities.data import AnnDataField, densify


class ModelInferenceService:
    """Service for model inference."""

    @staticmethod
    def _get_model_checkpoint_path(model_checkpoint_file_path: str) -> str:
        """
        Get model checkpoint path from GCS. It expects a filepath that doesn't include bucket protocol prefix and
        bucket name.

        :param model_checkpoint_file_path: Model checkpoint file path in GCS

        :return: Model checkpoint file path with GCS protocol prefix
        """
        return f"gs://{settings.PROJECT_BUCKET_NAME}/{model_checkpoint_file_path}"

    @classmethod
    @cache
    def _get_model_checkpoint_file(cls, model_file_path: str) -> t.BinaryIO:
        """
        Get model checkpoint from either local or GCS and load it using CellariumModule.

        :param model_file_path: Model checkpoint file path (from model db object)

        :return: CellariumModule object
        """
        model_checkpoint_path = cls._get_model_checkpoint_path(model_file_path)

        with open(model_checkpoint_path, "rb") as model_checkpoint_file:
            return BytesIO(model_checkpoint_file.read())

    @classmethod
    @cache
    def _load_module_from_checkpoint(cls, model_file_path: str) -> CellariumModule:
        """
        Load CellariumModule from checkpoint file.

        :param model_file_path: Model checkpoint file path (from model db object)

        :return: CellariumModule object
        """
        checkpoint_file = cls._get_model_checkpoint_file(model_file_path)

        return CellariumModule.load_from_checkpoint(checkpoint_file, map_location="cpu")

    @staticmethod
    def get_cache_info() -> dict[str, tuple[int, int, int | None, int]]:
        """
        Returns the cache info for the file and module cache.

        :return: A dict containing two entries: file_cache_info and module_cache_info.  Each
            contains a tuple with 4 values: hits, misses, maxsize, and currsize.
        """
        file_cache_info = ModelInferenceService._get_model_checkpoint_file.cache_info()
        module_cache_info = ModelInferenceService._load_module_from_checkpoint.cache_info()

        return {
            "file_cache_info": file_cache_info,
            "module_cache_info": module_cache_info,
        }

    @staticmethod
    def _create_cellarium_data_module(adata: anndata.AnnData) -> CellariumAnnDataDataModule:
        """
        Create CellariumAnnDataDataModule from anndata object that is ready to be used for model inference.

        :param adata: Anndata object

        :return: CellariumAnnDataDataModule
        """
        data_module = CellariumAnnDataDataModule(
            dadc=adata,
            batch_keys={
                "x_ng": AnnDataField(attr="X", convert_fn=densify),
                "var_names_g": AnnDataField(attr="var_names"),
                "total_mrna_umis_n": AnnDataField(attr="obs", key="total_mrna_umis"),
            },
            batch_size=len(adata),
            shuffle=False,
        )
        data_module.setup(stage="predict")
        return data_module

    def _get_output_from_model(self, model: models.CASModel, adata: anndata.AnnData) -> tuple[np.ndarray, list[str]]:
        """
        Get output from cellarium-ml model that predicts embeddings given an input adata.

        :param model: CAS Backend model db object
        :param adata: Object of :class:`anndata.AnnData` to embed.

        :return: Tuple of embeddings and obs_ids.
        """
        cellarium_module = ModelInferenceService._load_module_from_checkpoint(model.model_file_path)

        cellarium_data_module = self._create_cellarium_data_module(adata=adata)
        batch = next(iter(cellarium_data_module.predict_dataloader()))

        cellarium_output_dict = cellarium_module(batch)

        embeddings = cellarium_output_dict["x_ng"].numpy()

        obs_ids = adata.obs.index.tolist()
        return embeddings, obs_ids

    @staticmethod
    def _validate_model_output(embeddings: np.ndarray, obs_ids: list[str], model_info: models.CASModel) -> None:
        """
        Validate model output.

        :param embeddings: Embeddings
        :param obs_ids: List of observation ids
        :param model_info: CAS Backend model db object

        :raises ModelOutputError: If the length of obs_ids and embeddings are not the same or if the number of embedding
            dimensions is not equal to the model's embedding dimensions.
        """
        if embeddings.shape[0] != len(obs_ids):
            raise exceptions.ModelOutputError(
                f"The number of embeddings generated ({embeddings.shape[0]}) does not match "
                f"the number of observation IDs provided ({len(obs_ids)})."
            )

        if embeddings.shape[1] != model_info.embedding_dimension:
            raise exceptions.ModelOutputError(
                f"The dimensionality of the embeddings generated ({embeddings.shape[1]}) does not match "
                f"the expected embedding dimension ({model_info.embedding_dimension}) specified in model_info. "
                f"Ensure that the model is configured to produce embeddings of the correct dimensionality."
            )

    def embed_adata(self, adata: anndata.AnnData, model: models.CASModel) -> tuple[np.ndarray, list[str]]:
        """
        Embed adata using a specific model using Cellarium-ML model and pytorch.

        :param adata: Object of :class:`anndata.AnnData` to embed.
        :param model: Model object that contains relevant information to use for obtaining embedding.

        :return: ModelEmbeddings schema object.
        """
        embeddings, obs_ids = self._get_output_from_model(model=model, adata=adata)
        self._validate_model_output(embeddings=embeddings, obs_ids=obs_ids, model_info=model)

        return obs_ids, embeddings
