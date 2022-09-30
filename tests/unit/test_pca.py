"""
Tests of Incremental PCA implementation. We need to compare our implementation
with scikit-learn
"""

import unittest

import numpy as np
import torch

from sklearn import decomposition as sklearn_pca
from casp.ml.models import pca as casp_pca


class TestIncrementalPCA(unittest.TestCase):
    @staticmethod
    def create_test_data():
        np_batch1 = np.random.random((7, 4000)) * 10
        np_batch2 = np.random.random((7, 4000)) * 10
        np_batch3 = np.random.random((7, 4000)) * 10
        return np_batch1, np_batch2, np_batch3

    def test_casp_inc_pca_vs_sklearn_inc_pca(self):
        np_batch_1, np_batch_2, np_batch_3 = self.create_test_data()
        torch_batch_1 = torch.Tensor(np_batch_1)
        torch_batch_2 = torch.Tensor(np_batch_2)
        torch_batch_3 = torch.Tensor(np_batch_3)

        sk_learn_incremental_pca = sklearn_pca.IncrementalPCA(n_components=3)
        casp_incremental_pca = casp_pca.IncrementalPCA(n_components=3)

        sk_learn_incremental_pca.partial_fit(np_batch_1)
        sk_learn_incremental_pca.partial_fit(np_batch_2)
        sk_learn_incremental_pca.partial_fit(np_batch_3)

        casp_incremental_pca(torch_batch_1)
        casp_incremental_pca(torch_batch_2)
        casp_incremental_pca(torch_batch_3)

        singular_values_same = torch.all(
            torch.isclose(
                input=torch.Tensor(sk_learn_incremental_pca.singular_values_),
                other=casp_incremental_pca.singular_values,
                atol=1e-03,
            )
        ).item()

        components_same = torch.all(
            torch.isclose(
                input=torch.Tensor(sk_learn_incremental_pca.components_),
                other=casp_incremental_pca.components,
                atol=1e-03,
            )
        ).item()
        sk_transformed = sk_learn_incremental_pca.transform(np_batch_1)
        casp_transformed = casp_incremental_pca.transform(torch_batch_1)

        transformation_same = torch.all(
            torch.isclose(input=torch.Tensor(sk_transformed), other=casp_transformed, atol=1e-03)
        ).item()

        self.assertTrue(singular_values_same, msg="Singular Values are not the same as in scikit-learn")
        self.assertTrue(components_same, msg="PCA Components are not the same as in scikit-learn")
        self.assertTrue(transformation_same, msg="Transformation of input batch is not the same as in scikit-learn")


def main():
    """
    Standard entry point that forwards to unit test main.
    """
    unittest.main()


if __name__ == "__main__":
    main()
