import os
import random
import typing as t
import warnings
from collections import deque
from pathlib import Path

import anndata
import torch
from torch.utils.data import Dataset

from casp.ml import running_stats
from casp.ml.data import transforms

if t.TYPE_CHECKING:
    from google.cloud import storage


class CASDataset(Dataset):
    """
    Get anndata from GCP storage and convert it in torch tensors.
    """

    def _randomize_chunk_names(self) -> None:
        self._chunk_names_all = self._processed_chunks.copy()
        self._processed_chunks = []
        random.shuffle(self._chunk_names_all)

    def _get_random_chunk_name(self) -> str:
        try:
            name = self._chunk_names_all.pop()
        except IndexError:
            self._randomize_chunk_names()
            name = self._chunk_names_all.pop()
            self._epoch += 1
            if isinstance(self.transform, transforms.Compose):
                for transform in self.transform.transforms:
                    setattr(transform, "first_epoch", False)

            setattr(self.transform, "first_epoch", False)

        self._processed_chunks.append(name)
        return name

    @staticmethod
    def _get_local_anndata_path(filename) -> str:
        return f"./data/{filename.split('/')[-1]}"

    @property
    def bucket(self) -> "storage.Bucket":
        return self.storage_client.bucket(self.bucket_name)

    def _get_random_minibatch(self) -> anndata.AnnData:
        filename = self._get_random_chunk_name()
        blob = self.bucket.blob(filename)
        filepath = self._get_local_anndata_path(filename)
        blob.download_to_filename(filepath)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = anndata.read_h5ad(filepath)

        os.remove(filepath)
        return adata

    def _update_data(self) -> None:
        adata = self._get_random_minibatch()
        self.X = torch.Tensor(adata.raw.X.todense().astype(int))
        self.db_ids = torch.Tensor(adata.obs_names.values.astype(int))

        self.index = 0
        self.ids_q = deque(i for i in range(self.X.shape[0]))

    def __init__(
        self,
        bucket_name: str,
        storage_client: "storage.Client",
        use_gpu: bool = False,
        storage_path: str = "",
        chunk_size: int = 10000,
        transform: t.Optional["transforms.CASTransform"] = None,
        running_stat: t.Optional["running_stats.RunningStat"] = None,
    ) -> None:
        """
        :param storage_client: Google Storage client.
        :param bucket_name: A bucket which the data should be extracted
        :param use_gpu: True if data should be put to CUDA memory.
        :param chunk_size: Approximate size of the chunks
        :param transform: Transform applied to tensor
        """
        Path("./data").mkdir(parents=True, exist_ok=True)
        self.bucket_name = bucket_name
        self.chunk_size = chunk_size
        self.transform = transform
        self.running_stat = running_stat
        self.storage_client = storage_client
        self.storage_path = storage_path
        self.use_gpu = use_gpu
        self._epoch = 1
        self._chunk_names_all = [x.name for x in self.bucket.list_blobs(prefix=self.storage_path)]
        avoid_names = [self.storage_path, self.storage_path + "/"]
        self._chunk_names_all = list(filter(lambda x: x not in avoid_names, self._chunk_names_all))
        random.shuffle(self._chunk_names_all)
        self._processed_chunks = []
        self.X = None
        self.db_ids = None
        self.index = None
        self.ids_q = None
        self.apply_transform = running_stat is None
        self._update_data()

    def __len__(self) -> int:
        return (len(self._chunk_names_all) + len(self._processed_chunks)) * self.chunk_size

    def __getitem__(self, _: int) -> t.Tuple[torch.Tensor, int]:
        """
        :return: Tuple of tensor with a cell gene counts and db index.
        """
        try:
            index = self.ids_q.popleft()
        except IndexError:
            self._update_data()
            index = self.ids_q.popleft()

        x_i, db_index = self.X[index], self.db_ids[index]

        if self.use_gpu:
            x_i = x_i.cuda()
            db_index = db_index.cuda()

        if self.transform is not None and self.apply_transform:
            x_i = self.transform(x_i)

        return x_i, db_index

    def calculate_running_stats(self) -> None:
        """
        Calculate running stats. You probably will need to call this method before
        start training loop
        """
        for i in range(len(self)):
            if self._epoch == 1:
                x_i, db_index = self[i]
                self.running_stat(x_i)
        self.apply_transform = True
        print("Running Statistics are calculated")
