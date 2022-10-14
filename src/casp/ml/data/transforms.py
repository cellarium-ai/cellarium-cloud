import typing as t
import torch
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
    def __init__(self, C: int):
        self.C = C

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        return torch.log(1 + self.C * (tensor / torch.sum(tensor)))


class ColumnWiseNormalization(CASTransform):
    """
    It is required to run this normalization on 1 epoch to calculate the running
    statistics. This transform is applied only after the first epoch
    """

    def __init__(self):
        self.first_epoch = True
        self.running_sums: t.Optional[torch.Tensor] = None
        self.running_sums_squared: t.Optional[torch.Tensor] = None
        self.n_rows = 0
        self.mu: t.Optional[torch.Tensor] = None
        self.variance: t.Optional[torch.Tensor] = None

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.first_epoch:
            if self.running_sums is None:
                self.running_sums = torch.zeros(tensor.shape, device=tensor.device)
                self.running_sums_squared = torch.zeros(tensor.shape, device=tensor.device)

            self.running_sums += tensor
            self.running_sums_squared += torch.square(tensor)
            self.n_rows += 1
            return tensor
        else:
            if self.mu is None and self.variance is None:
                self.mu = self.running_sums / self.n_rows
                self.variance = self.running_sums_squared / (self.n_rows - 1) - (
                    self.n_rows / (self.n_rows - 1)
                ) * torch.square(self.mu)
                self.std = torch.sqrt(self.variance)
            res = (tensor - self.mu) / self.std
            res[torch.isnan(res)] = 0
            return res
