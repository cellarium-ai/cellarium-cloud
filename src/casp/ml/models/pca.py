import typing as t

import torch

from casp.ml.models.shared.pca import IncrementalPCABase
from casp.ml.models.shared.utils import svd_flip


class FullRankIncrementalPCA(IncrementalPCABase):
    def decompose(self, batch: torch.Tensor) -> t.Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        U, S, Vh = torch.linalg.svd(batch, full_matrices=False)
        U, Vh = svd_flip(U, Vh, U_decision=False)
        return U, S, Vh
