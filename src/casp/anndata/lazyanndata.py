from typing import Optional

import numpy as np
import braceexpand
from anndata import AnnData

from casp.anndata import read_h5ad_gcs


class LazyAnnData:
    """Lazy AnnData"""

    def __init__(
        self,
        urls,
        shard_size: int,
        last_shard_size: Optional[int] = None,
    ):
        self.shard_names = expand_urls(urls)
        assert isinstance(self.shard_names[0], str)
        self._shard = {"adata": None, "name": None}
        self.shard_size = shard_size
        self.last_shard_size = last_shard_size

    def __getitem__(self, index) -> AnnData:
        """Returns a sliced view of the object."""
        shard_idx, cell_idx = unpack_index(index, self.shard_size)
        adata = self._get_adata(shard_idx)
        return adata[cell_idx]

    def __len__(self) -> int:
        if last_shard_size is None:
            last_shard_size = len(read_h5ad_gcs(self.shard_names[-1]))
        return shard_size * (len(self.shard_names) - 1) + last_shard_size

    def _get_adata(self, shard_idx) -> AnnData:
        shard_name = self.shard_names[shard_idx]

        if shard_name != self._shard["name"]:
            self._update_shard(shard_name)

        return self._shard["adata"]

    def _update_shard(self, shard_name) -> None:
        adata = read_h5ad_gcs(shard_name)
        if shard_name != self.shard_names[-1]:
            assert len(adata) == self.shard_size
        self._shard["name"] = shard_name
        self._shard["adata"] = adata


def unpack_index(index: int, shard_size: int):
    """If index is list-like then all cells have to be from the same shard."""
    if isinstance(index, int):
        return divmod(index, shard_size)
    elif isinstance(index, (list, tuple)):
        shard_idx, cell_idx = np.divmod(index, shard_size)
        assert (shard_idx == shard_idx[0]).all(), "All cells have to be from the same shard"
        return shard_idx[0], cell_idx
    else:
        raise ValueError


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
