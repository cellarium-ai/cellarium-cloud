import typing as t

import anndata
import numpy as np
from cellarium.ml import CellariumAnnDataDataModule, CellariumModule
from smart_open import open

from casp.services import settings
from casp.services.db import models
from casp.services.model_inference import exceptions


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

    def _get_model_checkpoint_file(self, model: models.CASModel) -> t.BinaryIO:
        """
        Get model checkpoint from either local or GCS and load it using CellariumModule.

        :param model: Cellarium Cloud model db object

        :return: CellariumModule object
        """
        model_checkpoint_path = self._get_model_checkpoint_path(model.model_file_path)

        return open(model_checkpoint_path, "rb")

    def _get_output_from_model(
        self, model: models.CASModel, adata_file: t.BinaryIO
    ) -> t.Tuple[np.ndarray, t.List[str]]:
        """
        Get output from cellarium-ml model that predicts embeddings given an input adata file.

        :param model: Cellarium Cloud model db object
        :param adata_file: File object of :class:`anndata.AnnData` object to embed.

        :return: Tuple of embeddings and obs_ids.
        """
        adata = anndata.read_h5ad(adata_file)

        cellarium_checkpoint_file = self._get_model_checkpoint_file(model)
        cellarium_module = CellariumModule.load_from_checkpoint(cellarium_checkpoint_file, map_location="cpu")
        cellarium_checkpoint_file.seek(0)

        cellarium_data_module = CellariumAnnDataDataModule.load_from_checkpoint(
            cellarium_checkpoint_file, dadc=adata, batch_size=adata.n_obs, num_workers=0
        )
        cellarium_data_module.setup(stage="predict")
        batch = next(iter(cellarium_data_module.predict_dataloader()))

        cellarium_output_dict = cellarium_module(batch)

        embeddings = cellarium_output_dict["x_ng"].numpy()

        obs_ids = adata.obs.index.tolist()
        return embeddings, obs_ids

    @staticmethod
    def _validate_model_output(embeddings: np.ndarray, obs_ids: t.List[str], model_info: models.CASModel) -> None:
        """
        Validate model output.

        :param embeddings: Embeddings
        :param obs_ids: List of observation ids
        :param model_info: Cellarium Cloud model db object

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

    def embed_adata_file(self, file_to_embed: t.BinaryIO, model: models.CASModel) -> t.Tuple[np.ndarray, t.List[str]]:
        """
        Embed adata file using a specific model using Cellarium-ML model and pytorch.

        :param file_to_embed: File object of :class:`anndata.AnnData` object to embed.
        :param model: Model object that contains relevant information to use for obtaining embedding.

        :return: ModelEmbeddings schema object.
        """
        embeddings, obs_ids = self._get_output_from_model(model=model, adata_file=file_to_embed)
        self._validate_model_output(embeddings=embeddings, obs_ids=obs_ids, model_info=model)

        return obs_ids, embeddings
