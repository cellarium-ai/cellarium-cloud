from functools import cache
from io import BytesIO
import typing as t

import anndata
import numpy as np
from smart_open import open
import torch

from cellarium.cas_backend.apps.model_inference import exceptions
from cellarium.cas_backend.core.config import settings
from cellarium.cas_backend.core.db import models
from cellarium.ml import CellariumAnnDataDataModule, CellariumModule
from cellarium.ml.utilities.data import AnnDataField, densify


class ModelInferenceService:
    """Service for model inference."""

    LOGITS_OUTPUT_KEY = "y_logits_nc"
    LABEL_ATTRIBUTE_NAMES = ("active_cl_names", "y_categories")

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
        Load CellariumModule from checkpoint file and place it in eval mode.

        Lightning's ``load_from_checkpoint`` does not call ``.eval()`` and
        ``CellariumModule.configure_model`` does not either, so the module
        defaults to ``training=True``.  Calling ``.eval()`` here ensures that
        dropout layers are disabled and BatchNorm uses its running statistics
        (rather than per-batch statistics) for every request served from the
        cached module.

        :param model_file_path: Model checkpoint file path (from model db object)

        :return: CellariumModule in eval mode
        """
        checkpoint_file = cls._get_model_checkpoint_file(model_file_path)

        module = CellariumModule.load_from_checkpoint(checkpoint_file, map_location="cpu")
        cls._patch_filter_backward_compat(module)
        module.eval()
        return module

    @staticmethod
    def _patch_filter_backward_compat(module: CellariumModule) -> None:
        """Backfill ``ordering`` and ``allow_missing`` on ``Filter`` objects loaded from pre-0.0.8 checkpoints.

        cellarium-ml 0.0.8 added both ``Filter.ordering`` (default ``True``) and
        ``Filter.allow_missing`` (default ``False``).  Checkpoints saved with 0.0.7 or earlier
        pickle a ``Filter`` whose ``__dict__`` has neither key.  Both are plain Python bools so
        they are absent from the state dict and are not restored by ``load_state_dict``.
        """
        try:
            from cellarium.ml.transforms.filter import Filter
        except ImportError:
            return
        for m in module.modules():
            if isinstance(m, Filter):
                if not hasattr(m, "ordering"):
                    m.ordering = True
                if not hasattr(m, "allow_missing"):
                    m.allow_missing = False

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
        # This dummy value is only used by scVI and ignored by other models
        adata.obs["_scvi_batch_index_n"] = 0
        data_module = CellariumAnnDataDataModule(
            dadc=adata,
            batch_keys={
                "x_ng": AnnDataField(attr="X", convert_fn=densify),
                "var_names_g": AnnDataField(attr="var_names"),
                "total_mrna_umis_n": AnnDataField(attr="obs", key="total_mrna_umis"),
                "batch_index_n": AnnDataField(attr="obs", key="_scvi_batch_index_n"),
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

        with torch.inference_mode():
            cellarium_output_dict = cellarium_module(batch)

        embeddings = cellarium_output_dict["x_ng"].detach().numpy()

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

    def embed_adata(self, adata: anndata.AnnData, model: models.CASModel) -> tuple[list[str], np.ndarray]:
        """
        Embed adata using a specific model using Cellarium-ML model and pytorch.

        :param adata: Object of :class:`anndata.AnnData` to embed.
        :param model: Model object that contains relevant information to use for obtaining embedding.

        :return: Tuple of observation ids and embeddings.
        """
        embeddings, obs_ids = self._get_output_from_model(model=model, adata=adata)
        self._validate_model_output(embeddings=embeddings, obs_ids=obs_ids, model_info=model)

        return obs_ids, embeddings

    @staticmethod
    def _logits_to_probabilities(logits_nc: torch.Tensor) -> np.ndarray:
        """
        Convert raw classification logits into a per-cell probability distribution over active class labels.

        This is a pre-propagation distribution: each row sums to 1 and no mass sits on ancestor terms in the
        cell type ontology. Ancestor propagation is performed separately by
        :class:`~cellarium.cas_backend.apps.compute.services.annotation.classification.ClassificationOntologyUpliftStrategy`.

        :param logits_nc: Raw, pre-propagation logits of shape ``(n_cells, n_active_cats)``.

        :return: Row-normalized probabilities of shape ``(n_cells, n_active_cats)``.
        """
        return torch.nn.functional.softmax(logits_nc, dim=1).detach().numpy()

    @classmethod
    def _extract_class_labels(cls, cellarium_model: object, n_classes: int) -> list[str]:
        """
        Resolve the ordered class labels for a classification model's logits.

        :param cellarium_model: The underlying cellarium-ml model (``CellariumModule.model``).
        :param n_classes: Expected number of class labels (the width of the logits tensor).

        :raises exceptions.ModelOutputError: If none of ``LABEL_ATTRIBUTE_NAMES`` is present on
            ``cellarium_model``, or if the resolved label count does not match ``n_classes``.

        :return: Class labels in the same order as the logits columns.
        """
        for attribute_name in cls.LABEL_ATTRIBUTE_NAMES:
            if hasattr(cellarium_model, attribute_name):
                labels = [str(label) for label in getattr(cellarium_model, attribute_name)]
                break
        else:
            raise exceptions.ModelOutputError(
                f"Classification model {type(cellarium_model).__name__} does not expose any of the expected "
                f"label attributes {cls.LABEL_ATTRIBUTE_NAMES}."
            )

        if len(labels) != n_classes:
            raise exceptions.ModelOutputError(
                f"The number of class labels ({len(labels)}) does not match the number of logits columns "
                f"({n_classes})."
            )

        return labels

    def _get_classification_output_from_model(
        self, model: models.CASModel, adata: anndata.AnnData
    ) -> tuple[np.ndarray, list[str], list[str]]:
        """
        Get class probabilities from a cellarium-ml classification model given an input adata.

        :param model: CAS Backend model db object
        :param adata: Object of :class:`anndata.AnnData` to classify.

        :raises exceptions.ModelOutputError: If the model output does not contain ``y_logits_nc``.

        :return: Tuple of probabilities, class labels, and obs_ids.
        """
        cellarium_module = ModelInferenceService._load_module_from_checkpoint(model.model_file_path)

        cellarium_data_module = self._create_cellarium_data_module(adata=adata)
        batch = next(iter(cellarium_data_module.predict_dataloader()))

        with torch.inference_mode():
            output = cellarium_module(batch)

        if self.LOGITS_OUTPUT_KEY not in output:
            raise exceptions.ModelOutputError(
                f"Classification model output does not contain '{self.LOGITS_OUTPUT_KEY}'. "
                f"Output keys were: {sorted(output.keys())}."
            )

        probabilities = self._logits_to_probabilities(output[self.LOGITS_OUTPUT_KEY])
        labels = self._extract_class_labels(cellarium_module.model, n_classes=probabilities.shape[1])
        obs_ids = adata.obs.index.tolist()

        return probabilities, labels, obs_ids

    @staticmethod
    def _validate_classification_output(probabilities: np.ndarray, obs_ids: list[str], labels: list[str]) -> None:
        """
        Validate classification model output.

        :param probabilities: Per-cell class probabilities
        :param obs_ids: List of observation ids
        :param labels: Class labels corresponding to the probability columns

        :raises exceptions.ModelOutputError: If the number of rows does not match ``obs_ids`` or the number of
            columns does not match ``labels``.
        """
        if probabilities.shape[0] != len(obs_ids):
            raise exceptions.ModelOutputError(
                f"The number of probability rows generated ({probabilities.shape[0]}) does not match "
                f"the number of observation IDs provided ({len(obs_ids)})."
            )

        if probabilities.shape[1] != len(labels):
            raise exceptions.ModelOutputError(
                f"The number of probability columns generated ({probabilities.shape[1]}) does not match "
                f"the number of class labels ({len(labels)})."
            )

    def predict_adata(self, adata: anndata.AnnData, model: models.CASModel) -> tuple[list[str], np.ndarray, list[str]]:
        """
        Predict class probabilities for adata using a specific classification model.

        :param adata: Object of :class:`anndata.AnnData` to classify.
        :param model: Model object that contains relevant information to use for obtaining class probabilities.

        :return: Tuple of observation ids, class probabilities, and class labels.
        """
        probabilities, labels, obs_ids = self._get_classification_output_from_model(model=model, adata=adata)
        self._validate_classification_output(probabilities=probabilities, obs_ids=obs_ids, labels=labels)

        return obs_ids, probabilities, labels
