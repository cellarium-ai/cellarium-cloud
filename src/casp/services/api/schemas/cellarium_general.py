import datetime

from pydantic import BaseModel, Field


class CASModel(BaseModel):
    model_name: str = Field(example="human-pca-001")
    schema_name: str = Field(example="refdata-gex-mm10-2020-A")
    is_default_model: bool = Field(example=False)
    embedding_dimension: int = Field(example=512)

    class Config:
        orm_mode = True


class FeatureSchemaInfo(BaseModel):
    schema_name: str = Field(example="refdata-gex-mm10-2020-A")


class ApplicationInfo(BaseModel):
    application_version: str
    default_feature_schema: str = Field(example="refdata-gex-GRCh38-2020-A")


class UserInfo(BaseModel):
    username: str
    email: str
    should_ask_for_feedback: bool = Field(example=True)


class UserQuota(BaseModel):
    user_id: int
    quota: int
    remaining_quota: int
    quota_reset_date: datetime.datetime
