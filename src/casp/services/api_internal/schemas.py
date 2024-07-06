import typing as t

from pydantic import BaseModel, Field

from casp.services import settings


def convert_to_optional(schema):
    return {k: t.Optional[v] for k, v in schema.__annotations__.items()}


class CASModelInCreate(BaseModel):
    model_name: str = Field(example="my_model")
    model_file_path: str = Field(example="gs://my_bucket/my_model")
    embedding_dimension: int = Field(example=512)
    bq_dataset_name: str = Field(example="my_dataset", default=settings.DEFAULT_MODEL_BQ_DATASET_NAME)
    schema_name: str = Field(example="my_schema_name", default=settings.DEFAULT_SCHEMA_NAME)
    admin_use_only: bool = Field(example=True, default=True)
    is_default_model: bool = Field(example=False, default=False)

    class Config:
        orm_mode = True


class CASModelOut(CASModelInCreate):
    id: int = Field(example=1)
    bq_dataset_name: str = Field(example="my_dataset")
    schema_name: str = Field(example="my_schema_name")

    class Config:
        orm_mode = True


class CASIndexOut(BaseModel):
    index_name: str = Field(example="my_index")
    num_neighbors: int = Field(example=10)
    endpoint_id: str = Field(example="projects/1111111/locations/us-somewhere/indexEndpoints/11123212")
    embedding_dimension: int = Field(example=256)
    deployed_index_id: t.Optional[str] = Field(example="my_deployed_index_id")

    class Config:
        orm_mode = True


class CASIndexInCreate(CASIndexOut):
    model_name: str = Field(example="my_model")


class CASIndexInUpdate(BaseModel):
    num_neighbors: t.Optional[int] = Field(example=10)
    endpoint_id: t.Optional[str] = Field(example="projects/1111111/locations/us-somewhere/indexEndpoints/11123212")
    embedding_dimension: t.Optional[int] = Field(example=256)
    deployed_index_id: t.Optional[str] = Field(example="my_deployed_index_id")

    class Config:
        orm_mode = True
