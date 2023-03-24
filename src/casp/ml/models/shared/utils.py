import typing as t

import torch


def svd_flip(U, VT, U_decision=True) -> t.Tuple[torch.Tensor, torch.Tensor]:
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
