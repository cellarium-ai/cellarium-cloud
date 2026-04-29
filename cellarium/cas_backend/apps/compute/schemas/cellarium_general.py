import datetime

from pydantic import AliasPath, BaseModel, ConfigDict, Field, field_validator


class OntologicalColumnInfo(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    ontology_resource_name: str
    column_name: str
    description: str | None = None


class CASModel(BaseModel):
    model_config = ConfigDict(from_attributes=True, protected_namespaces=())

    model_name: str = Field(example="human-pca-001")
    description: str = Field(example="1st version of PCA model for human data", default="")
    schema_name: str = Field(example="refdata-gex-mm10-2020-A")
    is_default_model: bool = Field(example=False)
    embedding_dimension: int = Field(example=512)
    ontological_columns: list[OntologicalColumnInfo] = Field(
        default_factory=list,
        validation_alias=AliasPath("cell_info_metadata", "ontological_columns"),
    )

    # Set the default value for description if it is provided as None
    @field_validator("description", mode="before")
    @classmethod
    def set_description_default(cls, v):
        return v or ""


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
    weekly_quota: int
    remaining_weekly_quota: int
    quota_reset_date: datetime.datetime
    lifetime_quota: int | None
    remaining_lifetime_quota: int | None
    quota_increased: bool


class ClientVersionInput(BaseModel):
    client_version: str


class ClientVersionOutput(BaseModel):
    is_valid: bool
    min_version: str


class CellOntologyResourceResponse(BaseModel):
    cl_names: list[str]
    cell_ontology_term_id_to_cell_type: dict[str, str]
    children_dictionary: dict[str, list[str]]
    shortest_path_lengths_from_cell_root: dict[str, int]
    longest_path_lengths_from_cell_root: dict[str, int]
