import typing as t

from fastapi import APIRouter, Depends, File, Form, UploadFile

from casp.services.api import dependencies, schemas, services
from casp.services.db import models

cell_operations_router = APIRouter(prefix="/cellarium-cell-operations")

cell_operations_service = services.CellOperationsService()
cellarium_general_service = services.CellariumGeneralService()


@cell_operations_router.post(
    "/annotate", response_model=schemas.QueryAnnotationCellTypeSummaryStatisticsType, deprecated=True
)
async def annotate(
    file: UploadFile = File(),
    model_name: str = Form(),
    include_dev_metadata: t.Optional[bool] = Form(default=False),
    request_user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Deprecated endpoint. Will be removed in the future version.

    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param include_dev_metadata: Boolean flag indicating whether to include dev metadata in the response.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.

    :return: JSON response with annotations.
    """
    return await cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        user=request_user,
        file=file.file,
        model_name=model_name,
        include_extended_output=include_dev_metadata,
    )


@cell_operations_router.post(
    "/annotate-cell-type-summary-statistics-strategy",
    response_model=schemas.QueryAnnotationCellTypeSummaryStatisticsType,
)
async def annotate_cell_type_summary_statistics(
    file: UploadFile = File(),
    model_name: str = Form(),
    request_user: models.User = Depends(dependencies.authenticate_user),
    include_extended_output: t.Optional[bool] = Form(default=False),
):
    """
    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :param include_extended_output: Boolean flag indicating whether to include a breakdown by dataset id in the
        response.

    :return: JSON response with annotations.
    """
    return await cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        user=request_user,
        file=file.file,
        model_name=model_name,
        include_extended_output=include_extended_output,
    )


@cell_operations_router.post(
    "/annotate-cell-type-ontology-aware-strategy", response_model=schemas.QueryAnnotationOntologyAwareType
)
async def annotate_cell_type_ontology_aware_strategy(
    file: UploadFile = File(),
    model_name: str = Form(),
    prune_threshold: float = Form(),
    weighting_prefactor: float = Form(),
    request_user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Annotate a single anndata file with Cellarium CAS using the cell type statistics strategy. Input file should be
    validated and sanitized according to the model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :param prune_threshold: Prune threshold. Threshold for pruning the ontology graph in the response.
    :param weighting_prefactor: Distance exponential weighting prefactor.

    :return: JSON response with annotations.
    """

    return await cell_operations_service.annotate_cell_type_ontology_aware_strategy_with_activity_logging(
        user=request_user,
        file=file.file,
        model_name=model_name,
        prune_threshold=prune_threshold,
        weighting_prefactor=weighting_prefactor,
    )


@cell_operations_router.post(path="/nearest-neighbor-search", response_model=t.List[schemas.SearchQueryCellResult])
async def nearest_neighbor_search(
    file: UploadFile = File(),
    model_name: str = Form(),
    user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Search for similar cells in a single anndata file with Cellarium CAS. Input file should be validated and sanitized
    according to the model schema.
    Response represents a list of objects, each object contains query cell id and list of cas cell indices and distances
    that are nearest neighbors to the query cell.
    """
    return await cell_operations_service.search_adata_file(file=file.file, model_name=model_name, user=user)


@cell_operations_router.post(path="/query-cells-by-ids", response_model=t.List[schemas.CellariumCellMetadata])
def get_cells_by_ids(
    item: schemas.CellariumCellByIdsInput,
    user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Get cells by their ids from a single anndata file with Cellarium CAS. Input file should be validated and sanitized
    according to the model schema.
    """
    return cell_operations_service.get_cells_by_ids_for_user(
        user=user,
        cell_ids=item.cas_cell_ids,
        metadata_feature_names=item.metadata_feature_names,
    )
