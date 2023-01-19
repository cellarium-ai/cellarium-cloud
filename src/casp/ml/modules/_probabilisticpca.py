import math
from typing import Dict, Iterable, Optional, Sequence, Tuple, Union

import pyro
import pyro.distributions as dist
import torch
import torch.nn.functional as F
from pyro import poutine
from pyro.infer import Trace_ELBO
from pyro.nn import PyroModule, PyroParam

from scvi._constants import REGISTRY_KEYS
from scvi.module.base import PyroBaseModuleClass, auto_move_data
from scvi.nn import Encoder

from torch.distributions import constraints

_PROBABILISTIC_PCA_PYRO_MODULE_NAME = "probabilistic_pca"


class ProbabilisticPCAPyroModel(PyroModule):
    """
    A PyroModule that serves as the model for the ProbabilisticPCAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input features.
    n_components
        Number of components to model.
    """

    def __init__(self, n_obs: int, n_vars: int, n_components: int):
        super().__init__(_PROBABILISTIC_PCA_PYRO_MODULE_NAME)

        self.n_obs = n_obs
        self.n_vars = n_vars
        self.n_components = n_components
        # Populated by PyroTrainingPlan.

        self.W = PyroParam(lambda: torch.randn((n_components, n_vars)))
        self.sigma = PyroParam(lambda: torch.tensor(1.0), constraint=constraints.positive)

    @staticmethod
    def _get_fn_args_from_batch(tensor_dict: Dict[str, torch.Tensor]) -> Union[Iterable, dict]:

        x = tensor_dict[REGISTRY_KEYS.X_KEY]
        return (x,), {}

    @auto_move_data
    def forward(
        self,
        x: torch.Tensor
    ):
        """Forward pass."""
        with pyro.plate("cells", size=self.n_obs, subsample_size=x.shape[0]):
            z = pyro.sample(
                "z",
                dist.Normal(torch.zeros(self.n_components, device=x.device), 1).to_event(1),
            )
            pyro.sample(
                "feature_counts",
                dist.Normal(z @ self.W, self.sigma).to_event(1),
                obs=x,
            )


class ProbabilisticPCAPyroGuide(PyroModule):
    """
    A PyroModule that serves as the guide for the ProbabilisticPCAPyroModule class.

    Parameters
    ----------
    n_input
        Number of input features.
    n_components
        Number of components to model.
    """

    def __init__(self, n_obs: int, n_vars: int, n_components: int):
        super().__init__(_PROBABILISTIC_PCA_PYRO_MODULE_NAME)

        self.n_obs = n_obs
        self.n_vars = n_vars
        self.n_components = n_components

        self.L = PyroParam(lambda: torch.randn((n_vars, n_components)))
        self.z_sigma = PyroParam(lambda: torch.ones(n_components), constraint=constraints.positive)

    @auto_move_data
    def forward(
        self,
        x: torch.Tensor
    ):
        """Forward pass."""
        # Principal component distributions guide.
        with pyro.plate("cells", size=self.n_obs, subsample_size=x.shape[0]):
            z_mean = x @ self.L
            pyro.sample("z", dist.Normal(z_mean, self.z_sigma).to_event(1))


class ProbabilisticPCAPyroModule(PyroBaseModuleClass):
    """
    An amortized implementation of Latent Dirichlet Allocation :cite:p:`Blei03` implemented in Pyro.

    This module uses auto encoding variational Bayes to optimize the latent variables in the model.
    In particular, a fully-connected neural network is used as an encoder, which takes in feature counts
    as input and outputs the parameters of cell topic distribution. To employ the reparametrization trick
    stably, the Dirichlet priors are approximated by a Logistic-Normal distribution.
    The input feature counts tensor is a cell by features Bag-of-Words(BoW) representation
    of the counts. I.e. the model treats each cell's feature vector as ordered, not
    as unordered as in a Multinomial distribution.

    Parameters
    ----------
    n_input
        Number of input features.
    n_components
        Number of components to model.
    """

    def __init__(
        self,
        n_obs: int,
        n_vars: int,
        n_components: int,
    ):
        super().__init__()

        self.n_obs = n_obs
        self.n_vars = n_vars
        self.n_components = n_components

        self._model = ProbabilisticPCAPyroModel(self.n_obs, self.n_vars, self.n_components)
        self._guide = ProbabilisticPCAPyroGuide(self.n_obs, self.n_vars, self.n_components)
        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):  # noqa: D102
        return self._model

    @property
    def guide(self):  # noqa: D102
        return self._guide

    @auto_move_data
    @torch.inference_mode()
    def get_component_distribution(self, x: torch.Tensor) -> torch.Tensor:
        """
        Converts `x` to its inferred topic distribution.
        Parameters
        ----------
        x
            Counts tensor.
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.
        Returns
        -------
        A `x.shape[0] x n_topics` tensor containing the normalized topic distribution.
        """
        L = self._guide.L.detach()
        z_mean = x @ L
        return z_mean.cpu()

    @auto_move_data
    @torch.inference_mode()
    def get_elbo(self, x: torch.Tensor, library: torch.Tensor, n_obs: int) -> float:
        """
        Computes ELBO.

        Parameters
        ----------
        x
            Counts tensor.
        library
            Library sizes for each cell.
        n_obs
            Size of full batch. If n_obs < x.shape[0], ELBO is scaled by (n_obs / x.shape[0]).

        Returns
        -------
        The positive ELBO.
        """
        return Trace_ELBO().loss(self.model, self.guide, x, n_obs=n_obs)
