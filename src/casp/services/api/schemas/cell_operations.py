import typing as t

from pydantic import BaseModel, Field


class NeighborhoodCellTypeSummaryStatistics(BaseModel):
    cell_type: str = Field(example="erythrocyte")
    cell_count: int = Field(example=94)
    min_distance: float = Field(example=1589.847900390625)
    p25_distance: float = Field(example=1664.875244140625)
    median_distance: float = Field(example=1791.372802734375)
    p75_distance: float = Field(example=1801.3585205078125)
    max_distance: float = Field(example=1840.047119140625)


class CellTypeStatisticsExtendedObject(BaseModel):
    dataset_id: str = Field(example="a7a92fb49-50741b00a-244955d47")
    count_per_dataset: int = Field(example=10)
    min_distance: float = Field(example=1589.847900390625)
    max_distance: float = Field(example=1840.047119140625)
    median_distance: float = Field(example=1791.372802734375)
    mean_distance: float = Field(example=1791.372802734375)


class NeighborhoodCellTypeSummaryStatisticsExtended(NeighborhoodCellTypeSummaryStatistics):
    dataset_ids_with_counts: t.List[CellTypeStatisticsExtendedObject]


class QueryCellNeighborhoodAbstract(BaseModel):
    query_cell_id: str = Field(example="ATTACTTATTTAGTT-12311")
    matches: t.List


class QueryCellNeighborhoodCellTypeSummaryStatistics(QueryCellNeighborhoodAbstract):
    matches: t.List[NeighborhoodCellTypeSummaryStatistics]


class QueryCellNeighborhoodCellTypeSummaryStatisticsExtended(QueryCellNeighborhoodAbstract):
    query_cell_id: str = Field(example="ATTACTTATTTAGTT-12311")
    matches: t.List[NeighborhoodCellTypeSummaryStatisticsExtended]


class NeighborhoodSummaryOntologyAware(BaseModel):
    score: float = Field(example=0.789)
    cell_type_ontology_term_id: str = Field(example="CL:0000121")
    cell_type: str = Field(example="erythrocyte")


class QueryCellNeighborhoodOntologyAware(QueryCellNeighborhoodAbstract):
    query_cell_id: str = Field(example="ATTACTTATTTAGTT-12311")
    matches: t.List[NeighborhoodSummaryOntologyAware]
    total_weight: float = Field(example=11.23232)
    total_neighbors: int = Field(example=1023)
    total_neighbors_unrecognized: int = Field(example=5)


QueryAnnotationAbstractType: t.TypeAlias = t.List[QueryCellNeighborhoodAbstract]
QueryAnnotationOntologyAwareType: t.TypeAlias = t.List[QueryCellNeighborhoodOntologyAware]
QueryAnnotationCellTypeSummaryStatisticsType: t.TypeAlias = t.Union[
    t.List[QueryCellNeighborhoodCellTypeSummaryStatistics],
    t.List[QueryCellNeighborhoodCellTypeSummaryStatisticsExtended],
]
QueryAnnotationType: t.TypeAlias = t.Union[
    t.List[QueryCellNeighborhoodCellTypeSummaryStatisticsExtended],
    t.List[QueryCellNeighborhoodOntologyAware],
    t.List[QueryCellNeighborhoodCellTypeSummaryStatistics],
]


class SearchNeighborInfo(BaseModel):
    cas_cell_index: int = Field(example=123)
    distance: float = Field(example=0.123)


class SearchQueryCellResult(BaseModel):
    query_cell_id: str = Field(example="c")
    neighbors: t.List[SearchNeighborInfo]


class CellariumCellByIdsInput(BaseModel):
    cas_cell_ids: t.List[int]
    metadata_feature_names: t.List[str]


class CellariumCellMetadata(BaseModel):
    cas_cell_index: int = Field(default=None, example=123)
    cell_type: t.Optional[str] = Field(default=None, example="enterocyte")
    assay: t.Optional[str] = Field(default=None, example="10x 3' v2")
    disease: t.Optional[str] = Field(default=None, example="glioblastoma")
    donor_id: t.Optional[str] = Field(default=None, example="H20.33.013")
    is_primary_data: t.Optional[bool] = Field(default=None, example=True)
    development_stage: t.Optional[str] = Field(default=None, example="human adult stage")
    organism: t.Optional[str] = Field(default=None, example="Homo sapiens")
    self_reported_ethnicity: t.Optional[str] = Field(default=None, example="Japanese")
    sex: t.Optional[str] = Field(default=None, example="male")
    suspension_type: t.Optional[str] = Field(default=None, example="nucleus")
    tissue: t.Optional[str] = Field(default=None, example="cerebellum")
    total_mrna_umis: t.Optional[int] = Field(default=None, example=24312)

    # Ontology term IDs for the fields
    cell_type_ontology_term_id: t.Optional[str] = Field(default=None, example="CL:0000121")
    assay_ontology_term_id: t.Optional[str] = Field(default=None, example="EFO:0010550")
    disease_ontology_term_id: t.Optional[str] = Field(default=None, example="PATO:0000461")
    development_stage_ontology_term_id: t.Optional[str] = Field(default=None, example="HsapDv:0000053")
    organism_ontology_term_id: t.Optional[str] = Field(default=None, example="NCBITaxon:9606")
    self_reported_ethnicity_ontology_term_id: t.Optional[str] = Field(default=None, example="HANCESTRO:0019")
    sex_ontology_term_id: t.Optional[str] = Field(default=None, example="PATO:0000384")
    tissue_ontology_term_id: t.Optional[str] = Field(default=None, example="UBERON:0002037")
