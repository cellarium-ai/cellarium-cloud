"""
Regression tests for GitHub issue #230: eval mode in ModelInferenceService.

Both tests are designed to FAIL on the pre-fix code (module stays in training=True,
forward runs outside inference_mode) and PASS after the fix.
"""

from io import BytesIO
from unittest.mock import patch

import anndata
import pytest
import torch
import torch.nn as nn

from cellarium.cas_backend.apps.model_inference import services as model_inference_services
from cellarium.cas_backend.core.db import models
from cellarium.ml import CellariumModule
from tests.unit.fixtures import constants

# ---------------------------------------------------------------------------
# Stub module
# ---------------------------------------------------------------------------


class RecordingModule(nn.Module):
    """
    Minimal nn.Module stub whose forward records the mode state at call time.

    Using a real nn.Module subclass ensures that `.eval()` / `.train()` set
    the genuine ``self.training`` flag, and that ``__call__`` dispatches
    through the standard PyTorch hook chain into ``forward``.
    """

    def __init__(self, *, embedding_dim: int) -> None:
        super().__init__()
        self.embedding_dim = embedding_dim
        self.recorded_training: bool | None = None
        self.recorded_inference_mode: bool | None = None

    def forward(self, batch: dict) -> dict:
        """
        Record eval / inference-mode state, then return a plausible embedding tensor.

        :param batch: Batch dict produced by CellariumAnnDataDataModule.
        :return: Dict with key ``"x_ng"`` containing a zero tensor.
        """
        self.recorded_training = self.training
        self.recorded_inference_mode = torch.is_inference_mode_enabled()
        n_cells = batch["x_ng"].shape[0]
        return {"x_ng": torch.zeros(n_cells, self.embedding_dim)}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def cas_model() -> models.CASModel:
    """
    Minimal CASModel DB object for inference tests.

    Avoids pulling in the full ``populate_db`` fixture chain; only the fields
    read by ``_get_output_from_model`` and ``_validate_model_output`` are set.
    """
    return models.CASModel(
        id=1,
        model_name=constants.TEST_MODEL_NAME,
        admin_use_only=constants.TEST_MODEL_ADMIN_USE_ONLY,
        model_file_path=constants.TEST_MODEL_FILE_PATH,
        embedding_dimension=constants.TEST_EMBEDDING_DIMENSION,
    )


@pytest.fixture(autouse=True)
def clear_inference_caches():
    """
    Clear both classmethod caches before and after every test in this module.

    ``_load_module_from_checkpoint`` and ``_get_model_checkpoint_file`` are
    both decorated with ``@cache``, keyed on ``(cls, path)``.  Without this
    fixture a stub injected in one test leaks into the next, and a warm file
    cache hands out an exhausted BytesIO buffer to the second checkpoint load.
    """
    model_inference_services.ModelInferenceService._load_module_from_checkpoint.cache_clear()
    model_inference_services.ModelInferenceService._get_model_checkpoint_file.cache_clear()
    yield
    model_inference_services.ModelInferenceService._load_module_from_checkpoint.cache_clear()
    model_inference_services.ModelInferenceService._get_model_checkpoint_file.cache_clear()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_load_module_from_checkpoint_returns_eval_mode() -> None:
    """
    Module returned by ``_load_module_from_checkpoint`` must have ``training == False``.

    Patches at the level of ``CellariumModule.load_from_checkpoint`` so the
    real ``_load_module_from_checkpoint`` body runs — including the ``module.eval()``
    call introduced by the fix.  Old code (no ``.eval()``) → assertion fails.
    """
    stub = RecordingModule(embedding_dim=constants.TEST_EMBEDDING_DIMENSION)
    assert stub.training is True, "RecordingModule must start in training mode"

    with (
        patch.object(
            model_inference_services.ModelInferenceService,
            "_get_model_checkpoint_file",
            return_value=BytesIO(b""),
        ),
        patch.object(
            CellariumModule,
            "load_from_checkpoint",
            return_value=stub,
        ),
    ):
        module = model_inference_services.ModelInferenceService._load_module_from_checkpoint(
            constants.TEST_MODEL_FILE_PATH
        )

    assert module.training is False, (
        "_load_module_from_checkpoint must call .eval() before caching the module. "
        "Without it, dropout and BatchNorm run in training mode for every served request."
    )


def test_get_output_from_model_runs_in_eval_and_inference_mode(
    mock_valid_anndata: anndata.AnnData,
    cas_model: models.CASModel,
) -> None:
    """
    The forward pass inside ``_get_output_from_model`` must satisfy both:

    - ``module.training is False`` (eval mode, dropout/BN behave correctly)
    - ``torch.is_inference_mode_enabled() is True`` (no autograd, no BN stat mutation)

    Patches at the level of ``CellariumModule.load_from_checkpoint`` so the
    real ``_load_module_from_checkpoint`` body (including ``.eval()``) executes,
    and the real ``_create_cellarium_data_module`` and dataloader are exercised
    with ``mock_valid_anndata`` (20 cells × 500 genes, has ``total_mrna_umis``).
    """
    stub = RecordingModule(embedding_dim=constants.TEST_EMBEDDING_DIMENSION)

    with (
        patch.object(
            model_inference_services.ModelInferenceService,
            "_get_model_checkpoint_file",
            return_value=BytesIO(b""),
        ),
        patch.object(
            CellariumModule,
            "load_from_checkpoint",
            return_value=stub,
        ),
    ):
        service = model_inference_services.ModelInferenceService()
        embeddings, obs_ids = service._get_output_from_model(model=cas_model, adata=mock_valid_anndata)

    assert stub.recorded_training is False, (
        "Module must be in eval mode (training=False) during the forward pass. "
        "Old code left training=True, causing dropout to corrupt query embeddings."
    )
    assert stub.recorded_inference_mode is True, (
        "Forward pass must run inside torch.inference_mode(). "
        "Old code ran without it, allowing BatchNorm to mutate running stats across requests."
    )
    assert embeddings.shape == (20, constants.TEST_EMBEDDING_DIMENSION)
    assert obs_ids == [f"cell_{i}" for i in range(20)]
