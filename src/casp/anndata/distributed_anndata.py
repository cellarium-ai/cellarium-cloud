from collections import OrderedDict
from typing import List, Optional, Tuple, Union

import braceexpand
import numpy as np
import pandas as pd
from anndata import AnnData
from anndata._core.index import Index, _normalize_indices
from anndata.experimental.multi_files._anncollection import (AnnCollection,
                                                             AnnCollectionView,
                                                             ConvertType)

from casp.anndata import read_h5ad_gcs
from scvi.data._utils import get_anndata_attribute


class LRUCache:
    def __init__(self, maxsize: int):
        self.cache = OrderedDict()
        self.maxsize = maxsize

    def __getitem__(self, key):
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        raise KeyError(f"{key} not found in the cache")

    def __setitem__(self, key, value):
        while len(self.cache) >= self.maxsize:
            _, lru_value = self.cache.popitem(last=False)
            del lru_value
        self.cache[key] = value
        self.cache.move_to_end(key)

    def __contains__(self, key) -> bool:
        return key in self.cache


class LazyAnnData:
    """Lazy AnnData"""

    view_attrs = ["obs", "obsm", "layers", "var_names"]

    def __init__(self, filename, anncollection, cache, idx):
        self.filename = filename
        self.idx = idx
        self._n_obs = anncollection.shard_size
        self._obs_names = pd.RangeIndex(self.n_obs * self.idx, self.n_obs * (self.idx + 1))
        # reference anncollection
        self._anncollection = anncollection
        self.cache = cache

    @property
    def n_obs(self) -> int:
        return self._n_obs

    @property
    def n_vars(self) -> int:
        return len(self.var_names)

    @property
    def shape(self) -> Tuple[int, int]:
        return self.n_obs, self.n_vars

    @property
    def obs_names(self) -> pd.Index:
        return self._obs_names

    @property
    def cached(self) -> bool:
        return self.idx in self.cache

    @property
    def adata(self) -> AnnData:
        """Return backed AnnData from the filename"""
        if not self.cached:
            self.cache[self.idx] = read_h5ad_gcs(self.filename)
        return self.cache[self.idx]

    def __getattr__(self, attr):
        if attr in self.view_attrs:
            return getattr(self._anncollection.adatas[0].adata, attr)
        adata = self.adata
        if hasattr(adata, attr):
            return getattr(adata, attr)
        raise AttributeError(f"Backed AnnData object has no attribute '{attr}'")

    def __getitem__(self, idx) -> AnnData:
        adata = self.adata
        return adata[idx]

    def __repr__(self) -> str:
        if self.cached:
            buffered = "Buffered "
        else:
            buffered = ""
        backed_at = f" backed at {str(self.filename)!r}"
        descr = f"{buffered}LazyAnnData object with n_obs × n_vars = {self.n_obs} × {self.n_vars}{backed_at}"
        if self.cached:
            for attr in [
                "obs",
                "var",
                "uns",
                "obsm",
                "varm",
                "layers",
                "obsp",
                "varp",
            ]:
                keys = getattr(self, attr).keys()
                if len(keys) > 0:
                    descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr


class DistributedAnnCollection(AnnCollection):
    """Distributed AnnData Collection with Lazy Attributes"""

    is_view = False
    isbacked = False

    def __init__(
        self,
        filenames,
        shard_size: Optional[int] = None,
        last_shard_size: Optional[int] = None,
        indices_strict: bool = True,
        maxsize: int = 2,
        convert: Optional[ConvertType] = None,
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
        self.cache = LRUCache(maxsize=maxsize)
        self.maxsize = maxsize
        self.adatas = [LazyAnnData(filename, self, self.cache, idx) for idx, filename in enumerate(self.filenames)]
        self.uns = OrderedDict()
        super().__init__(
            adatas=self.adatas,
            join_obs=None,
            join_obsm=None,
            join_vars=None,
            label="dataset",
            keys=None,
            index_unique=None,
            convert=convert,
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
        assert len(indices) <= self.maxsize
        adatas = [None] * len(indices)
        for i, idx in enumerate(indices):
            if self.adatas[idx].cached:
                adatas[i] = self.adatas[idx].adata

        for i, idx in enumerate(indices):
            if not self.adatas[idx].cached:
                adatas[i] = self.adatas[idx].adata
        return adatas


@get_anndata_attribute.register
def _(
    adata: AnnCollection,
    attr_name: str,
    attr_key: Optional[str],
    mod_key: Optional[str] = None,
) -> Union[np.ndarray, pd.DataFrame]:
    return adata.lazy_attr(attr_name, attr_key)


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
