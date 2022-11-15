import typing as t

import torch

from casp.ml.utils import PickleMixin


class RunningStat(PickleMixin):
    def __call__(self, tensor):
        raise NotImplementedError

    def closure(self) -> None:
        return


class RunningSums(RunningStat):
    def __init__(self):
        self.running_sums: t.Optional[torch.Tensor] = None
        self.running_sums_squared: t.Optional[torch.Tensor] = None
        self.n = 0

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.running_sums is None or self.running_sums_squared is None:
            self.running_sums = torch.zeros(tensor.shape, device=tensor.device)
            self.running_sums_squared = torch.zeros(tensor.shape, device=tensor.device)

        self.running_sums += tensor
        self.running_sums_squared += torch.square(tensor)
        self.n += 1
        return tensor


class OnePassMeanVarStd(RunningSums):
    def __init__(self):
        super().__init__()
        self.mu: t.Optional[torch.Tensor] = None
        self.var: t.Optional[torch.Tensor] = None
        self.std: t.Optional[torch.Tensor] = None

    def closure(self) -> None:
        self.mu = self.running_sums / self.n
        self.var = self.running_sums_squared / (self.n - 1) - (self.n / (self.n - 1)) * torch.square(self.mu)
        self.std = torch.sqrt(self.var)
