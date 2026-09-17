"""
Regression tests for GitHub issue #230: eval mode in ModelInferenceService.

Both tests are designed to FAIL on the pre-fix code (module stays in training=True,
forward runs outside inference_mode) and PASS after the fix.
"""

from io import BytesIO
from unittest.mock import patch

import anndata
import numpy as np
import pytest
import torch
import torch.nn as nn

from cellarium.cas_backend.apps.model_inference import exceptions
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


class RecordingClassifierModule(nn.Module):
    """
    Minimal nn.Module stub for the classification path whose forward records the mode state at call time
    and returns a raw-logits classification output shaped like `SOCAM.predict`.

    Exposes `active_cl_names` and a `.model` property returning `self`, matching how `CellariumModule.model`
    resolves to the underlying cellarium-ml model in production. `model` is a property rather than
    `self.model = self` in `__init__`: `nn.Module.__setattr__` would otherwise register `self` as its own
    submodule, making `.modules()` / `.state_dict()` / `.to()` / `repr()` recurse infinitely over the
    self-referential module graph. Optionally includes an already-propagated `cell_type_probs_nc` key,
    mirroring SOCAM's real output shape, to prove the inference service never reads it.
    """

    def __init__(
        self,
        *,
        logits: torch.Tensor,
        active_cl_names: list[str],
        cell_type_probs: torch.Tensor | None = None,
    ) -> None:
        super().__init__()
        self.active_cl_names = active_cl_names
        self.recorded_training: bool | None = None
        self.recorded_inference_mode: bool | None = None
        self._logits = logits
        self._cell_type_probs = cell_type_probs

    @property
    def model(self) -> "RecordingClassifierModule":
        return self

    def forward(self, batch: dict) -> dict:
        """
        Record eval / inference-mode state, then return a plausible raw-logits classification output.

        :param batch: Batch dict produced by CellariumAnnDataDataModule.
        :return: Dict with key ``"y_logits_nc"``, and ``"cell_type_probs_nc"`` when configured.
        """
        self.recorded_training = self.training
        self.recorded_inference_mode = torch.is_inference_mode_enabled()
        output = {"y_logits_nc": self._logits}
        if self._cell_type_probs is not None:
            output["cell_type_probs_nc"] = self._cell_type_probs
        return output


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


def test_logits_to_probabilities_matches_classification_golden_probabilities() -> None:
    """
    `_logits_to_probabilities` softmaxes raw logits exactly as the classification golden fixture expects.

    Pins the same numbers from both sides: these logits are exactly what produces the `probabilities`
    array used by `tests/unit/test_classification_strategy.py`.
    """
    logits = torch.tensor([[3.0, 1.0, 0.0], [0.0, 3.0, 0.0]])
    probabilities = model_inference_services.ModelInferenceService._logits_to_probabilities(logits)

    expected = np.array(
        [
            [0.8437947344813395, 0.11419519938459449, 0.04201006613406605],
            [0.04527850074362907, 0.9094429985127419, 0.04527850074362907],
        ]
    )
    assert probabilities == pytest.approx(expected)


def test_extract_class_labels_from_active_cl_names() -> None:
    """`_extract_class_labels` reads `active_cl_names` when present (e.g. SOCAM)."""

    class FakeModel:
        active_cl_names = ["CL:0000789", "CL:0000084"]

    labels = model_inference_services.ModelInferenceService._extract_class_labels(FakeModel(), n_classes=2)
    assert labels == ["CL:0000789", "CL:0000084"]


def test_extract_class_labels_from_y_categories() -> None:
    """`_extract_class_labels` falls back to `y_categories` (e.g. LogisticRegression), coercing to str."""

    class FakeModel:
        y_categories = np.array(["CL:0000789", "CL:0000084"])

    labels = model_inference_services.ModelInferenceService._extract_class_labels(FakeModel(), n_classes=2)
    assert labels == ["CL:0000789", "CL:0000084"]


def test_extract_class_labels_raises_when_no_label_attribute() -> None:
    """`_extract_class_labels` raises `ModelOutputError` when neither label attribute is present."""

    class FakeModel:
        pass

    with pytest.raises(exceptions.ModelOutputError):
        model_inference_services.ModelInferenceService._extract_class_labels(FakeModel(), n_classes=2)


def test_extract_class_labels_raises_on_count_mismatch() -> None:
    """`_extract_class_labels` raises `ModelOutputError` when the label count disagrees with `n_classes`."""

    class FakeModel:
        active_cl_names = ["CL:0000789", "CL:0000084"]

    with pytest.raises(exceptions.ModelOutputError):
        model_inference_services.ModelInferenceService._extract_class_labels(FakeModel(), n_classes=3)


def test_get_classification_output_from_model_runs_in_eval_and_inference_mode(
    mock_valid_anndata: anndata.AnnData,
    cas_model: models.CASModel,
) -> None:
    """
    The forward pass inside `_get_classification_output_from_model` must satisfy both eval mode and
    inference mode -- the classification-path analogue of
    `test_get_output_from_model_runs_in_eval_and_inference_mode`, and the regression guard that keeps
    issue #230's fix from being bypassed by the classification path.
    """
    logits = torch.zeros(20, 2)
    stub = RecordingClassifierModule(logits=logits, active_cl_names=["CL:0000789", "CL:0000084"])

    with (
        patch.object(
            model_inference_services.ModelInferenceService,
            "_get_model_checkpoint_file",
            return_value=BytesIO(b""),
        ),
        patch.object(CellariumModule, "load_from_checkpoint", return_value=stub),
    ):
        service = model_inference_services.ModelInferenceService()
        probabilities, labels, obs_ids = service._get_classification_output_from_model(
            model=cas_model, adata=mock_valid_anndata
        )

    assert stub.recorded_training is False
    assert stub.recorded_inference_mode is True
    assert probabilities.shape == (20, 2)
    assert labels == ["CL:0000789", "CL:0000084"]
    assert obs_ids == [f"cell_{i}" for i in range(20)]


def test_get_classification_output_from_model_ignores_cell_type_probs_nc(
    mock_valid_anndata: anndata.AnnData,
    cas_model: models.CASModel,
) -> None:
    """
    When a checkpoint returns both `y_logits_nc` and an already-propagated `cell_type_probs_nc` (the real
    SOCAM output shape), `_get_classification_output_from_model` must use only the raw logits, never the
    pre-propagated probabilities -- using the latter would double-count ancestor propagation.
    """
    logits = torch.zeros(20, 3)
    logits[:, 0] = 3.0
    deliberately_different_probs = torch.full((20, 3), 1.0 / 3.0)
    stub = RecordingClassifierModule(
        logits=logits,
        active_cl_names=["CL:0000789", "CL:0000084", "CL:9999999"],
        cell_type_probs=deliberately_different_probs,
    )

    with (
        patch.object(
            model_inference_services.ModelInferenceService,
            "_get_model_checkpoint_file",
            return_value=BytesIO(b""),
        ),
        patch.object(CellariumModule, "load_from_checkpoint", return_value=stub),
    ):
        service = model_inference_services.ModelInferenceService()
        probabilities, _labels, _obs_ids = service._get_classification_output_from_model(
            model=cas_model, adata=mock_valid_anndata
        )

    expected = torch.nn.functional.softmax(logits, dim=1).numpy()
    assert probabilities == pytest.approx(expected)
    assert not np.allclose(probabilities, deliberately_different_probs.numpy())


def test_get_classification_output_from_model_raises_for_embedding_checkpoint(
    mock_valid_anndata: anndata.AnnData,
    cas_model: models.CASModel,
) -> None:
    """
    Registering an embedding checkpoint as `classification` must fail loudly rather than mis-serve:
    `_get_classification_output_from_model` raises `ModelOutputError` mentioning the missing `y_logits_nc`
    key.
    """
    stub = RecordingModule(embedding_dim=constants.TEST_EMBEDDING_DIMENSION)

    with (
        patch.object(
            model_inference_services.ModelInferenceService,
            "_get_model_checkpoint_file",
            return_value=BytesIO(b""),
        ),
        patch.object(CellariumModule, "load_from_checkpoint", return_value=stub),
    ):
        service = model_inference_services.ModelInferenceService()
        with pytest.raises(exceptions.ModelOutputError, match="y_logits_nc"):
            service._get_classification_output_from_model(model=cas_model, adata=mock_valid_anndata)


def test_patch_filter_backward_compat_backfills_ordering_and_allow_missing() -> None:
    """
    ``_patch_filter_backward_compat`` must backfill both ``ordering=True`` and
    ``allow_missing=False`` on any ``Filter`` submodule missing those attributes
    (i.e. loaded from a pre-0.0.8 checkpoint where neither existed).
    """
    from cellarium.ml.transforms.filter import Filter

    stub = RecordingModule(embedding_dim=constants.TEST_EMBEDDING_DIMENSION)

    # Attach a Filter child without ordering to simulate a pre-0.0.8 checkpoint.
    old_filter = Filter.__new__(Filter)
    # nn.Module.__init__ must be called to set up internal dicts, but we skip
    # Filter.__init__ so self.ordering is never set — reproducing the pickle state.
    nn.Module.__init__(old_filter)
    old_filter.filter_list = np.array(["GENE1", "GENE2"])
    assert not hasattr(old_filter, "ordering"), "Pre-condition: ordering must be absent"

    stub.old_filter = old_filter  # registers as a submodule via nn.Module.__setattr__

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

    for m in module.modules():
        if isinstance(m, Filter):
            assert hasattr(m, "ordering"), "ordering must be backfilled"
            assert m.ordering is True, "backfilled ordering must default to True"
            assert hasattr(m, "allow_missing"), "allow_missing must be backfilled"
            assert m.allow_missing is False, "backfilled allow_missing must default to False"
