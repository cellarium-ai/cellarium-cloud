import typing as t

if t.TYPE_CHECKING:
    import torch
    from cellarium.ml.module.incremental_pca import IncrementalPCA


class SCVIDIncrementalPCAWrapper:
    def __init__(self, model: "IncrementalPCA", embedding_dimension: int):
        """
        Wrapper around `scvid.module.incremental_pca.IncrementalPCA` for compatibility with
        CAS model inference service.

        :param model: Torch model
        :param embedding_dimension: Dimension of embedded data
        """
        self.embedding_dimension = embedding_dimension
        self.model = model

    def transform(self, batch: "torch.Tensor") -> "torch.Tensor":
        embeddings = self.model.predict(batch)
        return embeddings[:, : self.embedding_dimension]


# %%
