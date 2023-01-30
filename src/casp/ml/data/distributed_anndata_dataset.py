import typing as t

from torch import Tensor
from torch.utils.data import Dataset

from . import DistributedAnnCollection


class DistributedAnnCollectionDataset(Dataset):
    def __init__(self, dac: DistributedAnnCollection) -> None:
        self.dac = dac

    def __len__(self) -> int:
        return len(self.dac)

    def __getitem__(self, index: int) -> t.Tuple[Tensor, int]:
        """
        :return: Tuple of tensor with a cell gene counts and db index.
        """

        # TODO: review with yerdos
        v = self.dac[index]
        x_i = Tensor(v.X.todense().astype(int))
        db_index = Tensor(v.obs_names.values.astype(int))

        return x_i, db_index
