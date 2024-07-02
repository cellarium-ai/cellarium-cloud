import datetime

from pydantic import BaseModel, Field, validator


class CASModel(BaseModel):
    model_name: str = Field(example="human-pca-001")
    description: str = Field(example="1st version of PCA model for human data", default="")
    schema_name: str = Field(example="refdata-gex-mm10-2020-A")
    is_default_model: bool = Field(example=False)
    embedding_dimension: int = Field(example=512)

    # Set the default value for description if it is provided as None
    @validator('description', pre=True, always=True)
    def set_description_default(cls, v):
        return v or ""

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
