import typing as t
from pydantic import BaseModel


class ModelInferenceOutputBase(BaseModel):
    sample_ids: t.List[str]


class RepresentationModelOutput(ModelInferenceOutputBase):
    embeddings: t.Iterable


class ClassificationModelOutput(ModelInferenceOutputBase):
    probabilities: t.Iterable
    labels: t.Sequence[str]
