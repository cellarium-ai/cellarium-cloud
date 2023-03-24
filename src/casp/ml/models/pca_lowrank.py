import typing as t
import torch
from casp.ml.models.shared.pca import IncrementalPCABase
from casp.ml.models.shared.utils import svd_flip


class LowRankIncrementalPCA(IncrementalPCABase):
    def __init__(self, n_components, q=10000):
        self.q = q
        super().__init__(n_components=n_components)

    def decompose(self, batch: torch.Tensor) -> t.Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        U, S, Vh = torch.svd_lowrank(batch, q=self.q)
        U, Vh = svd_flip(U, Vh.T, U_decision=False)
        return U, S, Vh
