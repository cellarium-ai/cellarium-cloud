import typing as t
from functools import cache
from io import BytesIO

import anndata
import numpy as np
from cellarium.ml import CellariumAnnDataDataModule, CellariumModule
from smart_open import open

from casp.services import settings
from casp.services.db import models
from casp.services.model_inference import exceptions, schemas


class ModelInferenceServiceServiceInterface:
    @classmethod
    def run_inference(cls, adata: anndata.AnnData, model: models.CASModel) -> schemas.ModelInferenceOutputBase:
        """
        Run inference using the provided model on the given AnnData object. This method must be implemented by all
        subclasses.

        :param adata: Input annotated data matrix to run inference on.
        :param model: Cellarium Cloud model database object containing model metadata and checkpoint path.

        :return: Inference results wrapped in a subclass of ``ModelInferenceOutputBase``.
        """
        raise NotImplementedError


class CheckpointLoaderMixin:
    """
    Mixin with helper functions to load the cellarium module and data module from GCS
    """

    @staticmethod
    def _get_model_checkpoint_path(model_checkpoint_file_path: str) -> str:
        """
        Get a model checkpoint path from GCS. It expects a filepath that doesn't include bucket protocol prefix and
        bucket name.

        :param model_checkpoint_file_path: Model checkpoint file path in GCS

        :return: Model checkpoint file path with GCS protocol prefix
        """
        return f"gs://{settings.PROJECT_BUCKET_NAME}/{model_checkpoint_file_path}"

    @classmethod
    @cache
    def _get_model_checkpoint_file(cls, model_file_path: str) -> t.BinaryIO:
        """
        Get a model checkpoint from either local or GCS and load it using CellariumModule.

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

    @classmethod
    def _load_data_module_from_checkpoint(
        cls, model_file_path: str, adata: anndata.AnnData
    ) -> CellariumAnnDataDataModule:
        cellarium_checkpoint_file = cls._get_model_checkpoint_file(model_file_path)
        cellarium_checkpoint_file.seek(0)
        cellarium_data_module = CellariumAnnDataDataModule.load_from_checkpoint(
            cellarium_checkpoint_file, dadc=adata, batch_size=adata.n_obs, num_workers=0, shuffle=False
        )
        cellarium_data_module.setup(stage="predict")
        return cellarium_data_module


class RepresentationModelInferenceService(CheckpointLoaderMixin, ModelInferenceServiceServiceInterface):
    """Service for representation model inference."""

    @classmethod
    def _get_output_from_model(cls, model: models.CASModel, adata: anndata.AnnData) -> t.Tuple[np.ndarray, t.List[str]]:
        """
        Get output from cellarium-ml model that predicts embeddings given an input adata.

        :param model: Cellarium Cloud model db object
        :param adata: Object of class:`anndata.AnnData` to embed.

        :return: Tuple of embeddings and sample_ids.
        """
        cellarium_module = cls._load_module_from_checkpoint(model.model_file_path)

        cellarium_data_module = cls._load_data_module_from_checkpoint(
            model_file_path=model.model_file_path, adata=adata
        )

        batch = next(iter(cellarium_data_module.predict_dataloader()))

        cellarium_output_dict = cellarium_module(batch)

        embeddings = cellarium_output_dict["x_ng"].numpy()

        sample_ids = adata.obs.index.tolist()
        return embeddings, sample_ids

    @staticmethod
    def _validate_model_output(embeddings: np.ndarray, sample_ids: t.List[str], model_info: models.CASModel) -> None:
        """
        Validate model output.

        :param embeddings: Embeddings
        :param sample_ids: List of observation ids
        :param model_info: Cellarium Cloud model db object

        :raises ModelOutputError: If the length of sample_ids and embeddings are different, or if the number of embedding
            dimensions is not equal to the model's embedding dimensions.
        """
        if embeddings.shape[0] != len(sample_ids):
            raise exceptions.ModelOutputError(
                f"The number of embeddings generated ({embeddings.shape[0]}) does not match "
                f"the number of observation IDs provided ({len(sample_ids)})."
            )

        if embeddings.shape[1] != model_info.embedding_dimension:
            raise exceptions.ModelOutputError(
                f"The dimensionality of the embeddings generated ({embeddings.shape[1]}) does not match "
                f"the expected embedding dimension ({model_info.embedding_dimension}) specified in model_info. "
                f"Ensure that the model is configured to produce embeddings of the correct dimensionality."
            )

    @classmethod
    def run_inference(cls, adata: anndata.AnnData, model: models.CASModel) -> schemas.RepresentationModelOutput:
        """
        Embed adata using a specific model using Cellarium-ML model and pytorch.

        :param adata: Object of class:`anndata.AnnData` to embed.
        :param model: Model object that contains relevant information to use for getting embedding.

        :return: ModelEmbeddings schema object.
        """
        embeddings, sample_ids = cls._get_output_from_model(model=model, adata=adata)
        cls._validate_model_output(embeddings=embeddings, sample_ids=sample_ids, model_info=model)
        return schemas.RepresentationModelOutput(embeddings=embeddings, sample_ids=sample_ids)


class ClassificationModelInferenceService(CheckpointLoaderMixin, ModelInferenceServiceServiceInterface):
    @classmethod
    @cache
    def _load_module_from_checkpoint(cls, model_file_path: str) -> CellariumModule:
        """
        Load CellariumModule from checkpoint file.

        :param model_file_path: Model checkpoint file path (from model db object)

        :return: CellariumModule object
        """
        checkpoint_file = cls._get_model_checkpoint_file(model_file_path)
        checkpoint = CellariumModule.load_from_checkpoint(checkpoint_file, map_location="cpu")
        return checkpoint

    @classmethod
    def _get_output_from_model(
        cls, model: models.CASModel, adata: anndata.AnnData
    ) -> t.Tuple[np.ndarray, t.List[str], t.List[str]]:
        """
        Get output from cellarium-ml model that predicts embeddings given an input adata.

        :param model: Cellarium Cloud model db object
        :param adata: Object of class:`anndata.AnnData` to classify.

        :return: Tuple of embeddings, categories and sample_ids.
        """

        if "cell_type_ontology_term_id" not in adata.obs:
            adata.obs["cell_type_ontology_term_id"] = "CL_DUMMY"

        cellarium_module = cls._load_module_from_checkpoint(model_file_path=model.model_file_path)

        cellarium_module.model.actual_categories = 670  # total cell type categories to predict
        # cellarium_module.model.actual_categories = 2914  # total cell type categories to predict
        # cellarium_module.model.valid_mask = read_pkl_from_gcs(
        #     "gs://cellarium-file-system-cas-archive/curriculum/lrexp_human_validation_split_20241126/shared_meta/2914_socam_mask.pkl"
        # )
        # cellarium_module.model.target_row_descendent_col_torch_tensor = read_pkl_from_gcs(
        #     "gs://cellarium-file-system-cas-archive/curriculum/lrexp_human_validation_split_20241126/shared_meta/2914_target_row_descendent_col_torch_tensor_lrexp_human.pkl"
        # )
        # cellarium_module.model.probability_propagation_flag = True  # used for both training and validation
        #
        # cellarium_module.model.y_categories = read_pkl_from_gcs(
        #     "gs://cellarium-file-system-cas-archive/curriculum/lrexp_human_validation_split_20241126/shared_meta/sorted_2914_cell_type_names.pkl"
        # )
        cellarium_data_module = cls._load_data_module_from_checkpoint(
            model_file_path=model.model_file_path, adata=adata
        )

        batch = next(iter(cellarium_data_module.predict_dataloader()))

        forward_output_dict = cellarium_module.forward(batch)

        probabilities = forward_output_dict["cell_type_probs_nc"].detach().numpy()

        labels = cellarium_module.model.y_categories.tolist()
        # labels = read_pkl_from_gcs(
        #     "gs://cellarium-file-system-cas-archive/curriculum/lrexp_human_validation_split_20241126/shared_meta/sorted_2914_cell_type_names.pkl"
        # ).tolist()

        sample_ids = adata.obs.index.tolist()
        return probabilities, labels, sample_ids

    @staticmethod
    def _validate_model_output(probabilities: np.ndarray, sample_ids: t.List[str]) -> None:
        """
        Validate model output.

        :param probabilities: Probabilities from a classifier model
        :param sample_ids: List of observation ids

        :raises ModelOutputError: If the length of sample_ids and probabilities are different
        """
        if probabilities.shape[0] != len(sample_ids):
            raise exceptions.ModelOutputError(
                "Model output probabilities dimension doesn't correspond to input data length"
            )

    @classmethod
    def run_inference(cls, adata: anndata.AnnData, model: models.CASModel) -> schemas.ClassificationModelOutput:
        """
        Classify adata using a specific model using classifier Cellarium-ML model and pytorch.

        :param adata: Object of class:`anndata.AnnData` to embed.
        :param model: Model object that contains relevant information to use for predicting classes.

        :return: ClassificationModelOutput schema object.
        """
        probabilities, labels, sample_ids = cls._get_output_from_model(model=model, adata=adata)
        cls._validate_model_output(probabilities=probabilities, sample_ids=sample_ids)

        return schemas.ClassificationModelOutput(probabilities=probabilities, labels=labels, sample_ids=sample_ids)
