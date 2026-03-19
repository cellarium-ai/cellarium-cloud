from pydantic import BaseModel, Field


class NeighborhoodCellTypeSummaryStatistics(BaseModel):
    cell_type: str = Field(example="erythrocyte")
    cell_count: int = Field(example=94)
    min_distance: float = Field(example=1589.847900390625)
    p25_distance: float = Field(example=1664.875244140625)
    median_distance: float = Field(example=1791.372802734375)
    p75_distance: float = Field(example=1801.3585205078125)
    max_distance: float = Field(example=1840.047119140625)


class QueryCellNeighborhoodAbstract(BaseModel):
    query_cell_id: str = Field(example="ATTACTTATTTAGTT-12311")
    matches: list


class QueryCellNeighborhoodCellTypeSummaryStatistics(QueryCellNeighborhoodAbstract):
    matches: list[NeighborhoodCellTypeSummaryStatistics]


class NeighborhoodSummaryOntologyAware(BaseModel):
    score: float = Field(example=0.789)
    cell_type_ontology_term_id: str = Field(example="CL:0000121")
    cell_type: str = Field(example="erythrocyte")


class QueryCellNeighborhoodOntologyAware(QueryCellNeighborhoodAbstract):
    query_cell_id: str = Field(example="ATTACTTATTTAGTT-12311")
    matches: list[NeighborhoodSummaryOntologyAware]
    total_weight: float = Field(example=11.23232)
    total_neighbors: int = Field(example=1023)
    total_neighbors_unrecognized: int = Field(example=5)


type QueryAnnotationAbstractType = list[QueryCellNeighborhoodAbstract]
type QueryAnnotationOntologyAwareType = list[QueryCellNeighborhoodOntologyAware]
type QueryAnnotationCellTypeSummaryStatisticsType = list[QueryCellNeighborhoodCellTypeSummaryStatistics]
type QueryAnnotationType = (
    list[QueryCellNeighborhoodOntologyAware] | list[QueryCellNeighborhoodCellTypeSummaryStatistics]
)


class SearchNeighborInfo(BaseModel):
    cas_cell_index: int = Field(example=123)
    distance: float = Field(example=0.123)


class SearchQueryCellResult(BaseModel):
    query_cell_id: str = Field(example="c")
    neighbors: list[SearchNeighborInfo]


class CellariumCellByIdsInput(BaseModel):
    cas_cell_ids: list[int]
    metadata_feature_names: list[str]


class CellariumCellMetadata(BaseModel):
    cas_cell_index: int = Field(default=None, example=123)
    cell_type: str | None = Field(default=None, example="enterocyte")
    assay: str | None = Field(default=None, example="10x 3' v2")
    disease: str | None = Field(default=None, example="glioblastoma")
    donor_id: str | None = Field(default=None, example="H20.33.013")
    is_primary_data: bool | None = Field(default=None, example=True)
    development_stage: str | None = Field(default=None, example="human adult stage")
    organism: str | None = Field(default=None, example="Homo sapiens")
    self_reported_ethnicity: str | None = Field(default=None, example="Japanese")
    sex: str | None = Field(default=None, example="male")
    suspension_type: str | None = Field(default=None, example="nucleus")
    tissue: str | None = Field(default=None, example="cerebellum")
    total_mrna_umis: int | None = Field(default=None, example=24312)

    # Ontology term IDs for the fields
    cell_type_ontology_term_id: str | None = Field(default=None, example="CL:0000121")
    assay_ontology_term_id: str | None = Field(default=None, example="EFO:0010550")
    disease_ontology_term_id: str | None = Field(default=None, example="PATO:0000461")
    development_stage_ontology_term_id: str | None = Field(default=None, example="HsapDv:0000053")
    organism_ontology_term_id: str | None = Field(default=None, example="NCBITaxon:9606")
    self_reported_ethnicity_ontology_term_id: str | None = Field(default=None, example="HANCESTRO:0019")
    sex_ontology_term_id: str | None = Field(default=None, example="PATO:0000384")
    tissue_ontology_term_id: str | None = Field(default=None, example="UBERON:0002037")


class CacheInfo(BaseModel):
    file_cache_info: tuple[int, int, int | None, int]
    module_cache_info: tuple[int, int, int | None, int]
