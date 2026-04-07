import typing as t

from fastapi import APIRouter, Depends, File, Form, UploadFile

from cellarium.cas_backend.apps.compute import dependencies, schemas, services
from cellarium.cas_backend.core.db import models

AuthUser = t.Annotated[models.User, Depends(dependencies.authenticate_user)]
CellOpsService = t.Annotated[services.CellOperationsService, Depends(dependencies.get_cell_operations_service)]

cell_operations_router = APIRouter(prefix="/cellarium-cell-operations")


@cell_operations_router.post(
    "/annotate", response_model=schemas.QueryAnnotationCellTypeSummaryStatisticsType, deprecated=True
)
async def annotate(
    request_user: AuthUser,
    cell_operations_service: CellOpsService,
    file: UploadFile = File(),
    model_name: str = Form(),
    include_dev_metadata: bool | None = Form(default=False),
):
    """
    Deprecated endpoint. Will be removed in the future version.

    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param include_dev_metadata: Boolean flag indicating whether to include dev metadata in the response.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :param cell_operations_service: Service controller with domain logic responsible for cell operations

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
    request_user: AuthUser,
    cell_operations_service: CellOpsService,
    file: UploadFile = File(),
    model_name: str = Form(),
):
    """
    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :param include_extended_output: Boolean flag indicating whether to include a breakdown by dataset id in the
        response.
    :param cell_operations_service: Service controller with domain logic responsible for cell operations

    :return: JSON response with annotations.
    """
    return await cell_operations_service.annotate_cell_type_summary_statistics_strategy_with_activity_logging(
        user=request_user,
        file=file.file,
        model_name=model_name,
    )


@cell_operations_router.post(
    "/annotate-cell-type-ontology-aware-strategy", response_model=schemas.QueryAnnotationOntologyAwareType
)
async def annotate_cell_type_ontology_aware_strategy(
    request_user: AuthUser,
    cell_operations_service: CellOpsService,
    file: UploadFile = File(),
    model_name: str = Form(),
    prune_threshold: float = Form(),
    weighting_prefactor: float = Form(),
    ontology_column_name: str = Form(default="cell_type"),
):
    """
    Annotate a single anndata file with Cellarium CAS using the cell type statistics strategy. Input file should be
    validated and sanitized according to the model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.
    :param prune_threshold: Prune threshold. Threshold for pruning the ontology graph in the response.
    :param weighting_prefactor: Distance exponential weighting prefactor.
    :param ontology_column_name: Name of the ontological column to use for annotation (e.g. ``cell_type``, ``disease``).
        Defaults to ``cell_type``.
    :param cell_operations_service: Service controller with domain logic responsible for cell operations

    :return: JSON response with annotations.
    """

    return await cell_operations_service.annotate_cell_type_ontology_aware_strategy_with_activity_logging(
        user=request_user,
        file=file.file,
        model_name=model_name,
        prune_threshold=prune_threshold,
        weighting_prefactor=weighting_prefactor,
        ontology_column_name=ontology_column_name,
    )


@cell_operations_router.get("/cache-info", response_model=schemas.CacheInfo)
def get_cache_info(
    user: AuthUser,
    cell_operations_service: CellOpsService,
):
    """
    Return the cache info for the model checkpoint file and module caches.

    :param user: Authorized user object obtained  by token from `Bearer` header.
    :param cell_operations_service: Service controller with domain logic responsible for cell operations
    """

    return cell_operations_service.get_cache_info(user=user)


@cell_operations_router.post(path="/nearest-neighbor-search", response_model=list[schemas.SearchQueryCellResult])
async def nearest_neighbor_search(
    user: AuthUser,
    cell_operations_service: CellOpsService,
    file: UploadFile = File(),
    model_name: str = Form(),
):
    """
    Search for similar cells in a single anndata file with Cellarium CAS. Input file should be validated and sanitized
    according to the model schema.
    Response represents a list of objects, each object contains query cell id and list of cas cell indices and distances
    that are nearest neighbors to the query cell.
    """
    return await cell_operations_service.search_adata_file(file=file.file, model_name=model_name, user=user)
