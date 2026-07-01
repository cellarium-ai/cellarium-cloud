"""
Unit tests for cellarium.cas_backend.scripts.create_vsindex (library) and
scripts/create_vsindex.py (CLI shim).

Covers pure helpers, retry behaviour, CLI argument handling, and the
pre-flight wipe logic in _create_index — all without a live TileDB
vector-search runtime.
"""

import io

from click.testing import CliRunner
import numpy as np
import pytest
import tenacity
import tiledb

from cellarium.cas_backend.scripts.create_vsindex import (
    _choose_partitions,
    _create_index,
    _normalize_in_place,
    _process_chunk,
    _read_csv_gz,
)
from cellarium.cas_backend.scripts.create_vsindex_cli import main

# Stubs — minimal replacements for unittest.mock objects


class _FakeCtxMgr:
    """Context-manager stub; returns io.BytesIO(content) on enter."""

    def __init__(self, content: bytes):
        self._content = content

    def __enter__(self):
        return io.BytesIO(self._content)

    def __exit__(self, *_):
        pass


def _make_open_stub(raw: bytes):
    """Return a callable that mimics open() but yields BytesIO(raw)."""

    def _open(path, mode):
        return _FakeCtxMgr(raw)

    return _open


class _FakeIndex:
    """Stub index whose update_batch pops from a list of effects and counts calls."""

    def __init__(self, effects):
        self._effects = list(effects)
        self.update_batch_call_count = 0

    def update_batch(self, **kwargs):
        self.update_batch_call_count += 1
        effect = self._effects.pop(0)
        if isinstance(effect, BaseException):
            raise effect


class _FakeIngest:
    """Stub for vs.ingest; counts calls and returns a throwaway sentinel."""

    def __init__(self):
        self.call_count = 0

    def __call__(self, **kwargs):
        self.call_count += 1
        return object()


class _FakeVFS:
    """Stub tiledb.VFS; records is_dir / remove_dir call arguments."""

    def __init__(self, exists: bool):
        self._exists = exists
        self.is_dir_args: list[str] = []
        self.remove_dir_args: list[str] = []

    def is_dir(self, path: str) -> bool:
        self.is_dir_args.append(path)
        return self._exists

    def remove_dir(self, path: str) -> None:
        self.remove_dir_args.append(path)


class _FakeCreateVsindex:
    """Stub for create_vsindex; captures kwargs on each call."""

    def __init__(self):
        self.call_count = 0
        self.called_kwargs: dict = {}

    def __call__(self, **kwargs):
        self.call_count += 1
        self.called_kwargs = kwargs


# _choose_partitions


def test_choose_partitions_minimum_is_4096():
    assert _choose_partitions(1) == 4096


def test_choose_partitions_rounds_to_nearest_power_of_two():
    # target = min(65536, 7_000_000 // 1000) = 7000
    # nearest powers: 4096 (dist=2904) and 8192 (dist=1192); 8192 is closer
    result = _choose_partitions(7_000_000)
    assert result & (result - 1) == 0, "result must be a power of two"
    assert result == 8192


def test_choose_partitions_clamps_to_max_partitions():
    result = _choose_partitions(1_000_000_000, max_partitions=4096)
    assert result <= 4096


# _normalize_in_place


def test_normalize_in_place_unit_vector_unchanged():
    arr = np.array([[1.0, 0.0]], dtype=np.float32)
    _normalize_in_place(arr)
    assert np.allclose(arr, [[1.0, 0.0]])


def test_normalize_in_place_zero_vector_no_nan():
    arr = np.zeros((1, 3), dtype=np.float32)
    _normalize_in_place(arr)
    assert not np.any(np.isnan(arr))
    assert not np.any(np.isinf(arr))


def test_normalize_in_place_produces_unit_norms():
    rng = np.random.default_rng(42)
    arr = rng.random((10, 4)).astype(np.float32)
    _normalize_in_place(arr)
    norms = np.linalg.norm(arr, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-5)


# _read_csv_gz


def test_read_csv_gz_returns_correct_shape(monkeypatch):
    """smart_open decompresses; the stub returns already-decompressed bytes."""
    import cellarium.cas_backend.scripts.create_vsindex as mod

    raw_csv = b"0,0.1,0.2\n1,0.3,0.4\n2,0.5,0.6\n"
    monkeypatch.setattr(mod, "open", _make_open_stub(raw_csv))

    df = _read_csv_gz("/fake/batch_0.csv.gz")

    assert df.shape == (3, 3)
    assert list(df.iloc[:, 0]) == [0, 1, 2]


# _process_chunk — retry behaviour


@pytest.fixture
def process_chunk_env(monkeypatch):
    """Patches shared by both _process_chunk tests."""
    import cellarium.cas_backend.scripts.create_vsindex as mod

    monkeypatch.setattr(_process_chunk.retry, "wait", tenacity.wait_none())
    monkeypatch.setattr(mod, "open", _make_open_stub(b"0,0.1,0.2\n1,0.3,0.4\n"))


def test_process_chunk_retries_on_transient_tiledb_error(process_chunk_env):
    """_process_chunk must retry once when update_batch raises TileDBError, then succeed."""
    index = _FakeIndex(
        effects=[
            tiledb.TileDBError("GCS: List objects failed — Retry policy exhausted"),
            None,
        ]
    )

    result = _process_chunk(
        index=index,
        chunk_indices=[0],
        embeddings_prefix="/fake",
        embedding_dim=2,
        normalize=False,
    )

    assert result == 2
    assert index.update_batch_call_count == 2


def test_process_chunk_does_not_retry_non_tiledb_error(process_chunk_env):
    """Non-TileDBError exceptions must propagate immediately without retry."""
    index = _FakeIndex(effects=[ValueError("unexpected")])

    with pytest.raises(ValueError):
        _process_chunk(
            index=index,
            chunk_indices=[0],
            embeddings_prefix="/fake",
            embedding_dim=2,
            normalize=False,
        )

    assert index.update_batch_call_count == 1


# _create_index — pre-flight wipe


_SENTINEL_METRIC = object()


@pytest.fixture
def create_index_patches(monkeypatch):
    """Patches and inputs shared by both _create_index tests."""
    import cellarium.cas_backend.scripts.create_vsindex as mod

    monkeypatch.setattr(mod, "_resolve_distance_metric", lambda _: _SENTINEL_METRIC)
    monkeypatch.setattr(_create_index.retry, "wait", tenacity.wait_none())
    fake_ingest = _FakeIngest()
    monkeypatch.setattr(mod.vs, "ingest", fake_ingest, raising=False)
    return {
        "mod": mod,
        "fake_ingest": fake_ingest,
        "ids": np.array([0, 1], dtype=np.uint64),
        "embeddings": np.zeros((2, 4), dtype=np.float32),
    }


def test_create_index_wipes_existing_dir_before_creation(create_index_patches, monkeypatch, tmp_path):
    """_create_index removes an existing dir and calls vs.ingest exactly once."""
    vfs = _FakeVFS(exists=True)
    monkeypatch.setattr(create_index_patches["mod"].tiledb, "VFS", lambda: vfs)
    index_path = str(tmp_path / "idx")

    _create_index(
        ids=create_index_patches["ids"],
        embeddings=create_index_patches["embeddings"],
        index_path=index_path,
        total_batches=1,
        training_batches_used=1,
        total_training_rows=2,
        distance_metric="COSINE",
        normalize=False,
        max_partitions=65536,
    )

    assert vfs.is_dir_args == [index_path]
    assert vfs.remove_dir_args == [index_path]
    assert create_index_patches["fake_ingest"].call_count == 1


def test_create_index_skips_wipe_when_dir_absent(create_index_patches, monkeypatch, tmp_path):
    """_create_index must not call remove_dir when is_dir returns False."""
    vfs = _FakeVFS(exists=False)
    monkeypatch.setattr(create_index_patches["mod"].tiledb, "VFS", lambda: vfs)
    index_path = str(tmp_path / "idx")

    _create_index(
        ids=create_index_patches["ids"],
        embeddings=create_index_patches["embeddings"],
        index_path=index_path,
        total_batches=1,
        training_batches_used=1,
        total_training_rows=2,
        distance_metric="COSINE",
        normalize=False,
        max_partitions=65536,
    )

    assert vfs.remove_dir_args == []
    assert create_index_patches["fake_ingest"].call_count == 1


# CLI — argument handling


def test_main_missing_required_args_exits_nonzero():
    result = CliRunner().invoke(main, [])
    assert result.exit_code != 0


def test_main_invalid_distance_metric_exits_nonzero(tmp_path):
    result = CliRunner().invoke(
        main,
        [
            "--embeddings-prefix",
            str(tmp_path / "emb"),
            "--index-path",
            str(tmp_path / "idx"),
            "--total-batches",
            "10",
            "--embedding-dim",
            "64",
            "--distance-metric",
            "INVALID",
        ],
    )
    assert result.exit_code != 0


def test_main_calls_create_vsindex_with_correct_kwargs(monkeypatch, tmp_path):
    """CLI correctly maps all flags to create_vsindex keyword arguments."""
    import cellarium.cas_backend.scripts.create_vsindex_cli as cli_mod

    captured = _FakeCreateVsindex()
    monkeypatch.setattr(cli_mod, "create_vsindex", captured)
    prefix = str(tmp_path / "emb")
    idx = str(tmp_path / "idx")

    result = CliRunner().invoke(
        main,
        [
            "--embeddings-prefix",
            prefix,
            "--index-path",
            idx,
            "--total-batches",
            "10",
            "--embedding-dim",
            "64",
            "--distance-metric",
            "L2",
            "--no-normalize",
            "--max-partitions",
            "8192",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured.call_count == 1
    assert captured.called_kwargs == {
        "embeddings_prefix": prefix,
        "index_path": idx,
        "total_batches": 10,
        "embedding_dim": 64,
        "training_sample_size": 5_000_000,
        "distance_metric": "L2",
        "normalize": False,
        "max_partitions": 8192,
        "update_chunk_size": 50,
    }
