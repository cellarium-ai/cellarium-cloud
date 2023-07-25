import typing as t

import torch

from casp.ml.running_stats import OnePassMeanVarStd
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
        return torch.log(1 + self.sc_rna_normalization_pseudocount * (tensor / tensor.sum(dim=-1, keepdim=True)))


class ColumnWiseNormalization(CASTransform):
    """
    It is required to run this normalization on 1 epoch to calculate the running
    statistics. This transform is applied only after the first epoch
    """

    def __init__(self, one_pass_mean_var_std: "OnePassMeanVarStd"):
        self.mean_var_std = one_pass_mean_var_std

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.mean_var_std.mu is None:
            self.mean_var_std.calculate()
        res = (tensor - self.mean_var_std.mu) / self.mean_var_std.std
        res[torch.isnan(res)] = 0
        return res


class DummyTransform(CASTransform):
    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        return tensor
