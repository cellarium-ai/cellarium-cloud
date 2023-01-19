import scvi
import scanpy as sc
import matplotlib.pyplot as plt
from casp.ml.data import DistributedAnnCollection, read_h5ad_gcs
from scipy.sparse import issparse
from anndata.experimental.multi_files import AnnCollection
import anndata
import numpy as np
from casp.ml.models import ProbabilisticPCA

adata = DistributedAnnCollection(
    # "gs://dsp-cell-annotation-service/benchmark_v1/benchmark_v1.{000..002}.h5ad",
    "data/hca_subsampled_20k.{000..002}.h5ad",
    maxsize=4,
    convert={
        "X": lambda a: a.toarray() if issparse(a) else a,  # densify .X
        # 'obs': lambda a: np.asarray(a, dtype='float32'), # change dtype for all keys of .obsm
        # 'obs': {"batch" : lambda c: c.astype(str)} # change type only for one key of .obs
    },
)

# indexing and lazy attributes
#  obs = adata.lazy_attr("obs", key="age_group")
#  adata[[0, 21000, 50000]].obs["age_group"]
#  x1 = adata[[0, 20200]]
#  x2 = adata[:3, ['AL627309.1', 'AL669831.2', 'AL669831.5']]

# calculate mean and std
mu = adata[:].X.mean(0)
std = adata[:].X.std(0)

# Probabilistic PCA
adata.convert = {"X": lambda a: (a.toarray() - mu) / std}
adata.cache.maxsize = 2
ProbabilisticPCA.setup_anndata(adata)

pca_model = ProbabilisticPCA(adata, n_components=100)

pca_model.train(batch_size=1000, max_epochs=1)
# pca_model.save("pca_model/", overwrite=True)

latent = pca_model.get_latent_representation(batch_size=1000)
breakpoint()
pass
