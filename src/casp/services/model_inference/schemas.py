import typing as t

from pydantic import BaseModel


class ModelInferenceOutputBase(BaseModel):
    sample_ids: t.List[str]


class RepresentationModelOutput(ModelInferenceOutputBase):
    embeddings: t.Any  # numpy array or array-like


class ClassificationModelOutput(ModelInferenceOutputBase):
    probabilities: t.Any  # numpy array or array-like
    labels: t.Sequence[str]
