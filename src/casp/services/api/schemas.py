import typing as t

from pydantic import BaseModel, Field


class MatchInfo(BaseModel):
    cell_type: str = Field(example="erythrocyte")
    cell_count: int = Field(example=94)
    min_distance: float = Field(example=1589.847900390625)
    p25_distance: float = Field(example=1664.875244140625)
    median_distance: float = Field(example=1791.372802734375)
    p75_distance: float = Field(example=1801.3585205078125)
    max_distance: float = Field(example=1840.047119140625)


class QueryCell(BaseModel):
    query_cell_id: str = Field(example="99543")
    matches: t.List[MatchInfo]


class FeatureSchemaInfo(BaseModel):
    schema_name: str = Field(example="refdata-gex-mm10-2020-A")


class ApplicationInfo(BaseModel):
    application_version: str
    default_feature_schema: str = Field(example="refdata-gex-GRCh38-2020-A")
