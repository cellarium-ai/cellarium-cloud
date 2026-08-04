"""
Test-harness conftest for scripts unit tests.

On platforms where tiledb-vector-search is unavailable (e.g. linux/aarch64), inject
minimal sys.modules stubs so create_vsindex can be imported without the native wheel.
This is a test-only shim; it does not affect the production module's import statement.

Specifically, thanks to `from __future__ import annotations`, all `vs.IVFFlatIndex`
annotations in create_vsindex.py are lazy strings and are never evaluated at import
time.  The only runtime tiledb reference during import is `tiledb.TileDBError` used by
the `_tiledb_retry` decorator — the stub below satisfies exactly that.
"""

import sys
import types


def _install_tiledb_stub() -> None:
    """Insert minimal tiledb stub when the real package is absent."""
    real = sys.modules.get("tiledb")
    if real is not None and hasattr(real, "TileDBError"):
        return  # real tiledb already present; nothing to do

    tdb = types.ModuleType("tiledb")
    tdb.TileDBError = type("TileDBError", (Exception,), {})  # type: ignore[attr-defined]
    tdb.VFS = None  # will be replaced by monkeypatch in individual tests  # type: ignore[attr-defined]

    vsmod = types.ModuleType("tiledb.vector_search")
    tdb.vector_search = vsmod  # type: ignore[attr-defined]

    sys.modules.setdefault("tiledb", tdb)
    sys.modules.setdefault("tiledb.vector_search", vsmod)


_install_tiledb_stub()
