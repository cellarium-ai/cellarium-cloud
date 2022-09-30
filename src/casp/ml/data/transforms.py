import typing as t
import torch


class CASTransform:
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

    def __init__(self, n_columns, n_rows):
        self.first_epoch = True
        self.n_columns = n_columns
        self.n_rows = n_rows
        self.running_sums = torch.zeros(n_columns)
        self.running_sum_squares = torch.zeros(n_columns)
        self.mu: t.Optional[torch.Tensor] = None
        self.variance: t.Optional[torch.Tensor] = None

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.first_epoch:
            self.running_sums += tensor
            self.running_sum_squares += torch.square(tensor)
        else:
            if self.mu is None and self.variance is None:
                self.mu = self.running_sums / self.n_rows
                self.variance = (1 / (self.n_rows - 1)) * (self.running_sum_squares - (self.n_rows * self.mu**2))
                self.std = torch.sqrt(self.variance)

            return (tensor - self.mu) / self.std

        return tensor
