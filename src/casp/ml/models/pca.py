import typing as t

import numpy as np
import torch

from casp.ml.utils import PickleMixin


def svd_flip(U, VT, U_decision=True):
    """
    Flips the signs of U and VT for SVD in order to force deterministic output.

    Follows Sklearn convention by looking at U's maximum in columns
    as default.
    """
    if U_decision:
        max_abs_cols = torch.argmax(torch.abs(U), dim=0)
        signs = torch.sign(U[max_abs_cols, torch.arange(U.shape[1])])
    else:
        # rows of v, columns of u
        max_abs_rows = torch.argmax(torch.abs(VT), dim=1)
        signs = torch.sign(VT[torch.arange(VT.shape[0]), max_abs_rows])

    U *= signs
    VT *= signs[:, None]
    return U, VT


class IncrementalPCA(PickleMixin):
    def __init__(self, n_components):
        self.batch_num = 1
        self.n_components = n_components
        self.n_samples_seen = 0
        self._running_sum: t.Optional[torch.Tensor] = None
        self._running_sum_sq: t.Optional[torch.Tensor] = None
        self._curr_mean: t.Optional[torch.Tensor] = None
        self.singular_values: t.Optional[torch.Tensor] = None
        self.components: t.Optional[torch] = None
        self.transforms = None

    def __call__(self, X):
        batch = torch.clone(X)
        if any((self._running_sum is None, self._running_sum_sq is None, self._curr_mean is None)):
            self._running_sum = torch.zeros(batch.shape[1], device=batch.device)
            self._running_sum_sq = torch.zeros(batch.shape[1], device=batch.device)
            self._curr_mean = torch.zeros(batch.shape[1], device=batch.device)

        batch_size = batch.shape[0]
        self._running_sum += torch.sum(batch, dim=0)
        n = self.n_samples_seen + batch_size
        self._running_sum_sq += torch.sum(torch.square(batch), dim=0)
        last_mean = self._curr_mean
        curr_mean = self._running_sum / n
        curr_var = self._running_sum_sq / (n - 1) - (n / (n - 1)) * torch.square(curr_mean)

        if self.n_samples_seen == 0:
            batch -= curr_mean
        else:
            batch_mean = torch.mean(batch, dim=0)
            batch -= batch_mean
            mean_correction = np.sqrt((self.n_samples_seen / n) * batch_size) * (last_mean - batch_mean)
            batch = torch.vstack(
                (
                    self.singular_values.reshape((-1, 1)) * self.components,
                    batch,
                    mean_correction,
                )
            )

        U, S, Vh = torch.linalg.svd(batch, full_matrices=False)
        U, Vh = svd_flip(U, Vh, U_decision=False)

        explained_variance = S**2 / (n - 1)
        explained_variance_ratio = S**2 / torch.sum(curr_var * n)

        self.components = Vh[: self.n_components]
        self.singular_values = S[: self.n_components]
        self.explained_variance = explained_variance[: self.n_components]
        self.explained_variance_ratio = explained_variance_ratio[: self.n_components]
        self._curr_mean = curr_mean
        self._curr_var = curr_var
        self.n_samples_seen = n

    def transform(self, batch: torch.Tensor) -> torch.Tensor:
        return torch.matmul(batch - self._curr_mean, self.components[: self.n_components].T)

    def save(self, file_name):
        if self.transforms is None:
            raise ValueError("You need to save the model along with the transform " "to use it later with the dataset")

        super().save(file_name=file_name)
