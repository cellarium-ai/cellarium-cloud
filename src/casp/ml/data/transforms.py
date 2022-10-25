import typing as t

import torch

from casp.ml.running_stats import RunningSums
from casp.ml.utils import PickleMixin


class CASTransform(PickleMixin):
    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError


class Compose(CASTransform):
    """
    Inspired by `torchvision.transforms.Compose`
    """

    def __init__(self, transforms: t.List[CASTransform]):
        self.transforms = transforms

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        for transform in self.transforms:
            tensor = transform(tensor)
        return tensor

    def __repr__(self) -> str:
        format_string = f"{self.__class__.__name__}("
        for transform in self.transforms:
            format_string += "\n"
            format_string += f"    {transform}"

        format_string += "\n)"
        return format_string


class RowWiseNormalization(CASTransform):
    def __init__(self, sc_rna_normalization_pseudocount: int):
        self.sc_rna_normalization_pseudocount = sc_rna_normalization_pseudocount

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        return torch.log(1 + self.sc_rna_normalization_pseudocount * (tensor / torch.sum(tensor)))


class ColumnWiseNormalization(CASTransform):
    """
    It is required to run this normalization on 1 epoch to calculate the running
    statistics. This transform is applied only after the first epoch
    """

    def __init__(self, running_sums: "RunningSums"):
        self.running_sums = running_sums
        self.mu: t.Optional[torch.Tensor] = None
        self.variance: t.Optional[torch.Tensor] = None

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.mu is None or self.variance is None:
            self.mu = self.running_sums.running_sums / self.running_sums.n
            self.variance = self.running_sums.running_sums_squared / (self.running_sums.n - 1) - (
                self.running_sums.n / (self.running_sums.n - 1)
            ) * torch.square(self.mu)
            self.std = torch.sqrt(self.variance)

        res = (tensor - self.mu) / self.std
        res[torch.isnan(res)] = 0
        return res
