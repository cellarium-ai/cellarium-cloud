import collections.abc
import logging
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField, NumericalVarField
from casp.ml.modules import ProbabilisticPCAPyroModule
from scvi.utils import setup_anndata_dsp

from scvi.model.base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class ProbabilisticPCA(PyroSviTrainMixin, BaseModelClass):
    """
    Amortized Latent Dirichlet Allocation :cite:p:`Blei03`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.AmortizedLDA.setup_anndata`.
    n_topics
        Number of topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_topic_prior
        Prior of cell topic distribution. If `None`, defaults to `1 / n_topics`.
    topic_feature_prior
        Prior of topic feature distribution. If `None`, defaults to `1 / n_topics`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.AmortizedLDA.setup_anndata(adata)
    >>> model = scvi.model.AmortizedLDA(adata)
    >>> model.train()
    >>> feature_by_topic = model.get_feature_by_topic()
    >>> adata.obsm["X_LDA"] = model.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        n_components: int = 100,
    ):
        # in case any other model was created before that shares the same parameter names.
        pyro.clear_param_store()

        super().__init__(adata)

        n_obs = self.summary_stats.n_cells
        n_vars = self.summary_stats.n_vars

        self.module = ProbabilisticPCAPyroModule(
            n_obs=n_obs,
            n_vars=n_vars,
            n_components=n_components,
        )

        # TODO check this
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        **kwargs,
    ) -> Optional[AnnData]:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False, correct_data_format=False),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Converts a count matrix to an inferred PCA components.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.

        Returns
        -------
        A `n_obs x n_topics` Pandas DataFrame containing the normalized estimate
        of the topic distribution for each observation.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        dl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        transformed_xs = []
        for tensors in dl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            transformed_xs.append(self.module.get_component_distribution(x))
        transformed_x = torch.cat(transformed_xs).numpy()

        return pd.DataFrame(
            data=transformed_x,
            index=adata.obs_names,
            columns=[f"component_{i}" for i in range(transformed_x.shape[1])],
        )

    def get_elbo(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> float:
        """
        Return the ELBO for the data.

        The ELBO is a lower bound on the log likelihood of the data used for optimization
        of VAEs. Note, this is not the negative ELBO, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        The positive ELBO.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        dl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        elbos = []
        for tensors in dl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            library = x.sum(dim=1)
            elbos.append(self.module.get_elbo(x, len(dl.indices)))
        return np.mean(elbos)
