from typing import Optional, Tuple

import numpy as np
import pandas as pd
import braceexpand
import torch
import scipy.sparse as sp_sparse
from anndata import AnnData
from anndata._core.aligned_mapping import Layers

from casp.anndata import read_h5ad_gcs


class DistributedX(sp_sparse.spmatrix):
    def __init__(self, parent):
        self.dataset = parent

    def __getitem__(self, idx):
        """Returns a sliced view of the object."""
        shard_idx, cell_idx = unpack_index(idx, self.dataset.shard_size)
        adata = self.dataset._get_adata(shard_idx)
        return adata.X[cell_idx]

    @property
    def X(self):
        return self.dataset._buffer["adata"].X

    def __getattr__(self, attr):
        if hasattr(self.X, attr):
            return getattr(self.X, attr)
        raise AttributeError

    def __repr__(self):
        descr = f"BufferedArray shard {self.dataset._buffer['filename']}\n"
        descr += repr(self.X)
        return descr

    @property
    def shape(self):
        return self.dataset.shape


class DistributedDataFrame:
    """Distributed Layers"""

    def __init__(self, parent, df):
        self._df = df
        self._parent = parent

    def __getattr__(self, attr):
        if hasattr(self._df, attr):
            return getattr(self._df, attr)
        raise AttributeError

    def __getitem__(self, key):
        return self._parent._buffer["adata"].obs[key]

    def __setitem__(self, key, value):
        raise NotImplementedError

    def __len__(self) -> int:
        return len(self._parent)


class DistributedAnnData:
    """Distributed AnnData"""

    def __init__(
        self,
        filenames,
        shard_size: Optional[int] = None,
        last_shard_size: Optional[int] = None,
    ):
        self.filenames = expand_urls(filenames)
        assert isinstance(self.filenames[0], str)
        self._buffer = {"adata": None, "filename": None, "idx": None}
        if shard_size is None:
            shard_size = len(self._get_adata(0))
        self.shard_size = shard_size
        self.last_shard_size = last_shard_size
        self.is_view = False

        if self.last_shard_size is None:
            self.last_shard_size = len(self._get_adata(-1))
        self._n_obs = self.shard_size * (len(self.filenames) - 1) + self.last_shard_size
        self._layers = Layers(self)
        self._obs = DistributedDataFrame(self, self._buffer["adata"].obs)

    def __getattr__(self, attr):
        if attr in [
            "uns",
            "n_vars",
            "var_names",
            "isbacked",
            "_gen_repr",
        ]:
            return getattr(self._buffer["adata"], attr)
        raise AttributeError

    def __repr__(self):
        return self._gen_repr(self.n_obs, self.n_vars)

    @property
    def layers(self):
        return self._layers

    @property
    def obs(self):
        return self._obs

    @property
    def X(self):
        return DistributedX(self)

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        return self._n_obs

    @property
    def shape(self) -> Tuple[int, int]:
        """Shape of data matrix (:attr:`n_obs`, :attr:`n_vars`)."""
        return self.n_obs, self.n_vars

    def __len__(self) -> int:
        return self._n_obs

    def __getitem__(self, index) -> AnnData:
        """Returns a sliced view of the object."""
        shard_idx, cell_idx = unpack_index(index, self.shard_size)
        adata = self._get_adata(shard_idx)
        return adata[cell_idx]

    def _get_adata(self, shard_idx) -> AnnData:
        assert shard_idx < len(self.filenames)
        filename = self.filenames[shard_idx]

        if filename != self._buffer["filename"]:
            self._update_buffer(shard_idx)

        return self._buffer["adata"]

    def _update_buffer(self, shard_idx) -> None:
        assert shard_idx < len(self.filenames)
        filename = self.filenames[shard_idx]
        print(f"DEBUG updating shard to {filename}")
        adata = read_h5ad_gcs(filename)
        if filename != self.filenames[-1] and hasattr(self, "shard_size"):
            assert len(adata) == self.shard_size
        self._buffer["adata"] = adata
        self._buffer["filename"] = filename
        self._buffer["idx"] = shard_idx


def unpack_index(index: int, shard_size: int):
    """If index is list-like then all cells have to be from the same shard."""
    if isinstance(index, int):
        return divmod(index, shard_size)
    # elif isinstance(index, (list, tuple)):
    else:
        shard_idx, cell_idx = np.divmod(index, shard_size)
        assert (shard_idx == shard_idx[0]).all(), "All cells have to be from the same shard"
        return shard_idx[0], cell_idx
    #  else:
    #      raise ValueError


# https://github.com/webdataset/webdataset/blob/ab8911ab3085949dce409646b96077e1c1448549/webdataset/shardlists.py#L25-L33
def expand_urls(urls):
    if isinstance(urls, str):
        urllist = urls.split("::")
        result = []
        for url in urllist:
            result.extend(braceexpand.braceexpand(url))
        return result
    else:
        return list(urls)
