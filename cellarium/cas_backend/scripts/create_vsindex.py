"""
Create a TileDB IVF_FLAT vector search index from .csv.gz embedding batches.

Operates in three phases:
  1. Training      — read batch files until training_sample_size vectors are collected,
                     then create the IVF_FLAT index from those vectors. With an allow-list,
                     only cells whose soma_joinid is in the list are kept, and the training
                     set is filled to exactly training_sample_size allowed cells.
  2. Streaming     — append all remaining batch files to the index in chunks, plus any allowed
                     cells left over from the boundary batch when an allow-list is used.
  3. Consolidation — merge updates, retrain centroids, vacuum.
                     Skipped if no vectors were written in the streaming phase.

Example::

    python scripts/create_vsindex.py \\
        --embeddings-prefix gs://bucket/path/transform \\
        --index-path gs://bucket/path/ivflat_vsindex.soma \\
        --total-batches 1234 \\
        --embedding-dim 64

Importable directly from an installed package::

    from cellarium.cas_backend.scripts.create_vsindex import create_vsindex
"""

from __future__ import annotations

import gc
import io
import logging
import time

from google.api_core import exceptions as gax_exceptions  # for GCS transport errors in reads
import numpy as np
import pandas as pd
import pyarrow.csv as pcsv
from smart_open import open  # noqa: A001 — shadows built-in; handles gs:// and auto-decompresses .gz
import tenacity
import tiledb  # for tiledb.TileDBError used in retry and consolidate
from tiledb import vector_search as vs

logging.basicConfig(level=logging.INFO)

# Retry for tiledb operations (update_batch, consolidate_updates, vacuum).
# Both production failures are tiledb.TileDBError wrapping GCS "Retry policy exhausted".
# Using TileDBError (not Exception) means deterministic failures — bad distance metric,
# a bug raising ValueError — surface immediately instead of burning 5 × 60 s.
# tenacity = "8.2.3" is already in [tool.poetry.dependencies].
_tiledb_retry = tenacity.retry(
    retry=tenacity.retry_if_exception_type(tiledb.TileDBError),
    wait=tenacity.wait_exponential(multiplier=1, min=10, max=60),
    stop=tenacity.stop_after_attempt(5),
    reraise=True,
    before_sleep=tenacity.before_sleep_log(logging.getLogger(__name__), logging.WARNING),
)

# Retry for GCS reads via smart_open (_read_csv_gz).
# smart_open raises OSError / IOError for transport failures; google-api-core
# exceptions may bubble through for GCS-specific errors.
# Scope is intentionally narrow: _read_csv_gz already falls back pyarrow→pandas
# for parse errors — those must NOT be retried.
# Exception types here are disjoint from _tiledb_retry (TileDBError), so there is
# no compounding when _read_csv_gz is called inside a _tiledb_retry-wrapped function.
_transport_retry = tenacity.retry(
    retry=tenacity.retry_if_exception_type(
        (OSError, IOError, gax_exceptions.ServerError, gax_exceptions.ServiceUnavailable)
    ),
    wait=tenacity.wait_exponential(multiplier=1, min=10, max=60),
    stop=tenacity.stop_after_attempt(5),
    reraise=True,
    before_sleep=tenacity.before_sleep_log(logging.getLogger(__name__), logging.WARNING),
)

# Pure helper functions (no logging, no side effects)


def _choose_partitions(n_total: int, max_partitions: int = 65536) -> int:
    """Pick IVF partition count ≈ n_total/1000, rounded to nearest power of two, clamped to [4096, max_partitions]."""
    target = max(4096, min(max_partitions, n_total // 1000))
    p = 1
    while p < target:
        p <<= 1
    return p if (p - target) <= (target - (p >> 1)) else (p >> 1)


@_transport_retry
def _read_csv_gz(path: str) -> pd.DataFrame:
    """Open a .csv.gz file (smart_open decompresses transparently) and return a DataFrame."""
    with open(path, "rb") as f:
        buf = io.BytesIO(f.read())
    try:
        return pcsv.read_csv(buf, pcsv.ReadOptions(autogenerate_column_names=True)).to_pandas()
    except Exception:
        buf.seek(0)
        return pd.read_csv(buf, header=None)


@_transport_retry
def _load_allowed_ids(path: str) -> set[int]:
    """
    Read the allow-list CSV (single 'soma_joinid' column, header included) into a set.

    Unlike _read_csv_gz (header-less embedding batches), this CSV has a header row so
    pandas infers it. smart_open transparently decompresses a .gz path. The set is built
    once here and reused for every per-batch membership check.
    """
    with open(path, "rb") as f:
        buf = io.BytesIO(f.read())
    df = pd.read_csv(buf)
    return set(df["soma_joinid"].to_numpy(dtype=np.int64).tolist())


def _normalize_in_place(arr: np.ndarray) -> None:
    """L2-normalize rows of arr in-place. Clamps norms to >= 1e-12 to avoid division by zero."""
    norms = np.linalg.norm(arr, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)
    arr /= norms


def _resolve_distance_metric(name: str) -> object:
    """
    Return the tiledb DistanceMetric enum value for *name* (case-insensitive).

    Import of _tiledbvspy is deferred here so the module loads cleanly on platforms
    where the native extension is unavailable (e.g. the aarch64 CI stub).
    ValueError surfaces immediately — not retried by _tiledb_retry.
    """
    from tiledb.vector_search import _tiledbvspy as vspy  # noqa: PLC0415 — intentional lazy import

    metrics: dict[str, object] = {
        "COSINE": vspy.DistanceMetric.COSINE,
        "L2": vspy.DistanceMetric.L2,
    }
    result = metrics.get(name.upper())
    if result is None:
        raise ValueError(f"Unsupported distance metric: {name!r}. Choose one of {list(metrics)}")
    return result


# Phase functions


def _load_training_data(
    embeddings_prefix: str,
    total_batches: int,
    embedding_dim: int,
    training_sample_size: int,
    normalize: bool,
    allowed_ids: set[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, int, int, tuple[np.ndarray, np.ndarray] | None]:
    """
    Phase 1 (loading): Read batch files until training_sample_size vectors are collected.

    When allowed_ids is supplied, each batch is filtered to allow-listed soma_joinids and the
    training set is filled to exactly training_sample_size allowed cells; the boundary batch's
    leftover allowed cells are returned as a carryover for the update set (never dropped).

    Returns (ids [uint64], embeddings [float32], training_batches_used, total_training_rows,
    carryover) where carryover is None or a (ids [uint64], embeddings [float32]) pair of raw
    (un-normalized) leftover cells.
    """
    t0 = time.time()
    logging.info("Phase 1: Loading training data …")

    train_dfs = []
    carry_df = None
    total_training_rows = 0
    training_batches_used = 0

    for i in range(total_batches):
        df = _read_csv_gz(f"{embeddings_prefix}/batch_{i}.csv.gz")
        training_batches_used = i + 1
        if allowed_ids is None:
            # No thinning: whole-batch fill — original behavior, unchanged.
            train_dfs.append(df)
            total_training_rows += len(df)
            if i % 100 == 0 and i > 0:
                logging.info(f"  Loaded {training_batches_used} batches, {total_training_rows} vectors so far …")
            if total_training_rows >= training_sample_size:
                break
        else:
            # Thinning: filter by allow-list, then fill to EXACTLY training_sample_size,
            # carrying the boundary batch's leftover allowed cells into the update set.
            df = df[df.iloc[:, 0].isin(allowed_ids)]
            if total_training_rows + len(df) >= training_sample_size:
                need = training_sample_size - total_training_rows
                train_dfs.append(df.iloc[:need])
                carry_df = df.iloc[need:]
                total_training_rows += need
                break
            train_dfs.append(df)
            total_training_rows += len(df)
            if i % 100 == 0 and i > 0:
                logging.info(f"  Loaded {training_batches_used} batches, {total_training_rows} vectors so far …")

    if total_training_rows == 0:
        raise ValueError(
            "No training vectors collected. When an allow-list is supplied, it must match at least "
            "one cell in the input batches."
        )

    ids = np.empty(total_training_rows, dtype=np.uint64)
    embeddings = np.empty((total_training_rows, embedding_dim), dtype=np.float32)

    offset = 0
    for df in train_dfs:
        n = len(df)
        ids[offset : offset + n] = df.iloc[:, 0].to_numpy(dtype=np.uint64)
        embeddings[offset : offset + n] = df.iloc[:, 1:].to_numpy(dtype=np.float32)
        offset += n

    carryover = None
    if carry_df is not None and len(carry_df):
        carry_ids = carry_df.iloc[:, 0].to_numpy(dtype=np.uint64)
        carry_emb = carry_df.iloc[:, 1:].to_numpy(dtype=np.float32)
        carryover = (carry_ids, carry_emb)

    del train_dfs
    gc.collect()

    if normalize:
        _normalize_in_place(embeddings)

    logging.info(
        f"Training data loaded: {total_training_rows} vectors from {training_batches_used} batches "
        f"in {time.time()-t0:.1f}s"
    )
    return ids, embeddings, training_batches_used, total_training_rows, carryover


@_tiledb_retry
def _create_index(
    ids: np.ndarray,
    embeddings: np.ndarray,
    index_path: str,
    total_batches: int,
    training_batches_used: int,
    total_training_rows: int,
    distance_metric: str,
    normalize: bool,
    max_partitions: int,
) -> vs.IVFFlatIndex:
    """
    Phase 1 (creation): Build the IVF_FLAT index from training vectors.

    Always wipes and rebuilds index_path from scratch (pre-flight wipe runs on every attempt,
    including retries, so partial state from a failed vs.ingest is cleared before each retry).

    Returns the created index object. Caller is responsible for freeing ids/embeddings.
    """
    t0 = time.time()

    # Pre-flight wipe — runs on EVERY attempt (initial call and each retry).
    _vfs = tiledb.VFS()
    if _vfs.is_dir(index_path):
        logging.warning(f"Removing existing data at {index_path} before index creation.")
        _vfs.remove_dir(index_path)

    estimated_total = int((total_training_rows / training_batches_used) * total_batches)
    partitions = _choose_partitions(estimated_total, max_partitions)

    dist_metric = _resolve_distance_metric(distance_metric)

    logging.info(f"Creating IVF_FLAT index: {partitions} partitions, ~{estimated_total} total vectors")

    index = vs.ingest(
        index_type="IVF_FLAT",
        index_uri=index_path,
        input_vectors=embeddings,
        external_ids=ids.astype(np.int64),
        partitions=partitions,
        training_sample_size=total_training_rows,
        distance_metric=dist_metric,
        normalized=normalize,
    )

    logging.info(f"Phase 1 complete: index created in {time.time()-t0:.1f}s")
    return index


def _pack_and_update(index: vs.IVFFlatIndex, ids: np.ndarray, emb: np.ndarray) -> None:
    """Pack row embeddings into a TileDB object array and write them via update_batch."""
    n = len(ids)
    emb_obj = np.empty(n, dtype=object)
    for j in range(n):
        emb_obj[j] = emb[j]
    index.update_batch(vectors=emb_obj, external_ids=ids.astype(np.int64))


@_tiledb_retry
def _process_chunk(
    index: vs.IVFFlatIndex,
    chunk_indices: list[int],
    embeddings_prefix: str,
    embedding_dim: int,
    normalize: bool,
    allowed_ids: set[int] | None = None,
) -> int:
    """
    Read a chunk of batches and write them to the index via update_batch.

    Returns chunk_rows (number of vectors written). Decorated with _tiledb_retry;
    retry-safe because the updates array allows_duplicates=False and a later timestamp
    supersedes the earlier write, so re-running the same chunk is idempotent.
    """
    chunk_dfs = [_read_csv_gz(f"{embeddings_prefix}/batch_{i}.csv.gz") for i in chunk_indices]
    if allowed_ids is not None:
        chunk_dfs = [df[df.iloc[:, 0].isin(allowed_ids)] for df in chunk_dfs]
    chunk_rows = sum(len(df) for df in chunk_dfs)
    if chunk_rows == 0:
        return 0

    chunk_ids = np.empty(chunk_rows, dtype=np.uint64)
    chunk_emb = np.empty((chunk_rows, embedding_dim), dtype=np.float32)

    offset = 0
    for df in chunk_dfs:
        n = len(df)
        chunk_ids[offset : offset + n] = df.iloc[:, 0].to_numpy(dtype=np.uint64)
        chunk_emb[offset : offset + n] = df.iloc[:, 1:].to_numpy(dtype=np.float32)
        offset += n

    del chunk_dfs
    gc.collect()

    if normalize:
        _normalize_in_place(chunk_emb)

    _pack_and_update(index, chunk_ids, chunk_emb)

    del chunk_ids, chunk_emb
    gc.collect()

    return chunk_rows


@_tiledb_retry
def _write_carryover(index: vs.IVFFlatIndex, ids: np.ndarray, emb: np.ndarray, normalize: bool) -> int:
    """
    Write the boundary batch's leftover allowed cells (carried from Phase 1) to the update set.

    Retry-safe for the same reason as _process_chunk (update_batch is idempotent under a later
    timestamp). Returns rows written.
    """
    if normalize:
        _normalize_in_place(emb)
    _pack_and_update(index, ids, emb)
    return len(ids)


def _stream_updates(
    index: vs.IVFFlatIndex,
    embeddings_prefix: str,
    embedding_dim: int,
    remaining_batches: list[int],
    normalize: bool,
    update_chunk_size: int,
    allowed_ids: set[int] | None = None,
    carryover: tuple[np.ndarray, np.ndarray] | None = None,
) -> int:
    """
    Phase 2: Stream remaining batches into the index in chunks.

    Writes any boundary-batch carryover first (leftover allowed cells from Phase 1), then the
    remaining batches (each filtered by allowed_ids when thinning). Returns the total number of
    vectors written across carryover and chunks.
    """
    total_updated = 0

    if carryover is not None:
        c_ids, c_emb = carryover
        if len(c_ids):
            logging.info(f"Phase 2: writing {len(c_ids)} carried-over boundary-batch cells …")
            total_updated += _write_carryover(index, c_ids, c_emb, normalize)

    if not remaining_batches:
        if total_updated == 0:
            logging.info("Phase 2: No remaining batches; skipping.")
        return total_updated

    t0 = time.time()
    logging.info(f"Phase 2: Streaming {len(remaining_batches)} update batches …")

    for chunk_num, chunk_start in enumerate(range(0, len(remaining_batches), update_chunk_size)):
        chunk_indices = remaining_batches[chunk_start : chunk_start + update_chunk_size]
        total_updated += _process_chunk(index, chunk_indices, embeddings_prefix, embedding_dim, normalize, allowed_ids)
        if chunk_num % 5 == 4:
            batches_done = chunk_start + len(chunk_indices)
            logging.info(
                f"  Phase 2: {batches_done}/{len(remaining_batches)} batches, {total_updated} vectors written …"
            )

    logging.info(f"Phase 2 complete: {total_updated} vectors in {time.time()-t0:.1f}s")
    return total_updated


@_tiledb_retry
def _consolidate(
    index: vs.IVFFlatIndex,
    index_path: str,
    normalize: bool,
) -> vs.IVFFlatIndex:
    """
    Phase 3: Consolidate updates, retrain centroids, vacuum, and reload the index.

    No pre-flight wipe is needed on retry. Unlike vs.ingest (which writes a new index
    group from scratch and may leave a corrupt partial group on failure), consolidate_updates
    works on the existing updates array: unconsolidated fragments remain there until vacuum
    removes them, so a retry simply reprocesses what is still pending — idempotent by design.

    Returns reloaded index with updated .size.
    """
    t0 = time.time()
    logging.info("Phase 3: Consolidating updates …")

    index.consolidate_updates(retrain_index=True, normalized=normalize)
    index.vacuum()

    index = vs.IVFFlatIndex(index_path)
    logging.info(f"Phase 3 complete: {index.size} vectors, consolidation in {time.time()-t0:.1f}s")
    return index


# Public orchestrator


def create_vsindex(
    embeddings_prefix: str,
    index_path: str,
    total_batches: int,
    embedding_dim: int,
    training_sample_size: int = 5_000_000,
    distance_metric: str = "COSINE",
    normalize: bool = True,
    max_partitions: int = 65536,
    update_chunk_size: int = 50,
    allow_list_csv: str | None = None,
) -> None:
    """
    Create a TileDB IVF_FLAT vector search index from .csv.gz embedding batches.
    Callable directly from pipeline code without going through the CLI.

    When allow_list_csv is given, only cells whose soma_joinid appears in that CSV are indexed;
    the training set is filled to exactly training_sample_size allowed cells and every remaining
    allowed cell (including the boundary batch's leftover) is streamed into the update set.
    """
    t0 = time.time()
    logging.info(
        f"Starting vsindex creation: {total_batches} batches, dim={embedding_dim}, "
        f"metric={distance_metric}, normalize={normalize}"
    )

    allowed_ids = _load_allowed_ids(allow_list_csv) if allow_list_csv else None
    if allowed_ids is not None:
        logging.info(f"Thinning enabled: {len(allowed_ids)} allowed cell IDs from {allow_list_csv}")

    ids, embeddings, training_batches_used, total_training_rows, carryover = _load_training_data(
        embeddings_prefix, total_batches, embedding_dim, training_sample_size, normalize, allowed_ids
    )
    index = _create_index(
        ids,
        embeddings,
        index_path,
        total_batches,
        training_batches_used,
        total_training_rows,
        distance_metric,
        normalize,
        max_partitions,
    )
    del ids, embeddings
    gc.collect()

    remaining_batches = list(range(training_batches_used, total_batches))
    total_updated = _stream_updates(
        index,
        embeddings_prefix,
        embedding_dim,
        remaining_batches,
        normalize,
        update_chunk_size,
        allowed_ids,
        carryover,
    )

    if total_updated > 0:
        index = _consolidate(index, index_path, normalize)

    logging.info(f"vsindex creation complete in {time.time()-t0:.1f}s — final size: {index.size}")
