from collections import OrderedDict
from typing import List, Optional

import braceexpand
import pandas as pd
from anndata import AnnData
from anndata._core.index import Index, _normalize_indices
from anndata.experimental.multi_files._anncollection import (AnnCollection,
                                                             AnnCollectionView)

from casp.anndata import read_h5ad_gcs


class LazyAnnData:
    """Lazy AnnData"""

    view_attrs = ["obs", "obsm", "layers", "var_names"]

    def __init__(self, filename, anncollection, idx):
        self.filename = filename
        self.idx = idx
        self._anncollection = anncollection
        self._n_obs = anncollection.shard_size
        self._buffer = anncollection._buffer
        self.buffer_size = anncollection.buffer_size
        self._obs_names = pd.RangeIndex(self.n_obs * self.idx, self.n_obs * (self.idx + 1))

    @property
    def n_obs(self):
        return self._n_obs

    @property
    def n_vars(self):
        return len(self.var_names)

    @property
    def shape(self):
        return self.n_obs, self.n_vars

    @property
    def obs_names(self):
        return self._obs_names

    def materialize(self):
        return self._anncollection.materialize(self.idx)[0]

    def __getattr__(self, attr):
        if attr in self.view_attrs:
            return getattr(self._anncollection.materialize(0)[0], attr)
        adata = self.materialize()
        if hasattr(adata, attr):
            return getattr(adata, attr)
        raise AttributeError

    def __getitem__(self, idx) -> AnnData:
        adata = self.materialize()
        return adata[idx]


class DistributedAnnCollection(AnnCollection):
    """Distributed AnnData Collection with Lazy Attributes"""

    def __init__(
        self,
        filenames,
        shard_size: Optional[int] = None,
        last_shard_size: Optional[int] = None,
        indices_strict: bool = True,
        buffer_size: int = 2,
    ):
        self.filenames = expand_urls(filenames)
        assert isinstance(self.filenames[0], str)
        first_adata = read_h5ad_gcs(self.filenames[0])
        if shard_size is None:
            shard_size = len(first_adata)
        self.shard_size = shard_size
        self.last_shard_size = last_shard_size
        if self.last_shard_size is None:
            self.last_shard_size = self.shard_size
            # self.last_shard_size = len(self._get_adata(-1))
        self._buffer = OrderedDict()
        self.buffer_size = buffer_size
        self.adatas = [
            LazyAnnData(filename, self, idx, first_adata) if idx == 0 else LazyAnnData(filename, self, idx)
            for idx, filename in enumerate(self.filenames)
        ]
        super().__init__(
            adatas=self.adatas,
            join_obs=None,
            join_obsm=None,
            join_vars=None,
            label="dataset",
            keys=None,
            index_unique=None,
            convert=None,
            harmonize_dtypes=False,
            indices_strict=True,
        )

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)
        adatas_indices = [i for i, e in enumerate(resolved_idx[0]) if e is not None]
        self.materialize(adatas_indices)

        return AnnCollectionView(self, self.convert, resolved_idx)

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"DistributedAnnCollection object with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        descr += f"\n  constructed from {len(self.filenames)} AnnData objects"
        for attr, keys in self._view_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    view of {attr}: {str(keys)[1:-1]}"
        for attr in self._attrs:
            keys = list(getattr(self, attr).keys())
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        if "obs" in self._view_attrs_keys:
            keys = list(self.obs.keys())
            if len(keys) > 0:
                descr += f"\n    own obs: {str(keys)[1:-1]}"

        return descr

    def materialize(self, indices) -> List[AnnData]:
        if isinstance(indices, int):
            indices = (indices,)
        assert len(indices) <= self.buffer_size
        adatas = [None] * len(indices)
        for i, idx in enumerate(indices):
            if idx in self._buffer:
                self._buffer.move_to_end(idx)
                adatas[i] = self._buffer[idx]

        for i, idx in enumerate(indices):
            if idx not in self._buffer:
                while len(self._buffer) >= self.buffer_size:
                    _, adata = self._buffer.popitem(last=False)
                    del adata
                self._buffer[idx] = read_h5ad_gcs(self.filenames[idx])
                adatas[i] = self._buffer[idx]
        return adatas


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
