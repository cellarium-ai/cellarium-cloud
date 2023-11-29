import typing as t

from pydantic import BaseModel, Field


class ModelEmbeddings(BaseModel):
    """Model embedding schema"""

    obs_ids: t.List[str] = Field(
        ...,
        example=["ACTCGAGATACGATCTC", "ACTGCGCTAGCTA", "ACTCGAGATCG"],
        description="Cell ids from the input file index",
    )
    embeddings_b64: str = Field(
        ...,
        example="AAAAAA==",
        description="Cell embeddings in base64 format. Could be decoded with `base64.b64decode` into numpy",
    )
