import typing as t

from fastapi import APIRouter, Depends, File, Form, UploadFile

from casp.services.api import dependencies, schemas, services
from casp.services.api.services.consensus_engine import ConsensusStrategyType
from casp.services.db import models

cell_operations_router = APIRouter(prefix="/cellarium-cell-operations")

cell_operations_service = services.CellOperationsService()
cellarium_general_service = services.CellariumGeneralService()


@cell_operations_router.post("/annotate", response_model=schemas.QueryAnnotationCellTypeCountType)
async def annotate(
    file: UploadFile = File(),
    model_name: str = Form(),
    include_dev_metadata: t.Optional[bool] = Form(default=None),
    request_user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param include_dev_metadata: Boolean flag indicating whether to include dev metadata in the response.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.

    :return: JSON response with annotations.
    """

    return await cell_operations_service.annotate_adata_file(
        user=request_user,
        file=file.file,
        model_name=model_name,
        include_dev_metadata=include_dev_metadata,
        consensus_strategy=ConsensusStrategyType.CELL_TYPE_COUNT,
    )


@cell_operations_router.post(
    "/annotate-ontology-aware-strategy",
    response_model=schemas.QueryAnnotationOntologyAwareType,
)
async def annotate_ontology_aware_strategy(
    file: UploadFile = File(),
    model_name: str = Form(),
    normalize: t.Optional[bool] = Form(default=None),
    request_user: models.User = Depends(dependencies.authenticate_user),
):
    """
    Annotate a single anndata file with Cellarium CAS. Input file should be validated and sanitized according to the
    model schema.

    :param file: Byte object of :class:`anndata.AnnData` file to annotate.
    :param model_name: Model name to use for annotation. See `/list-models` endpoint for available models.
    :param normalize: Boolean flag indicating whether to normalize the data before annotation.
    :param request_user: Authorized user object obtained  by token from `Bearer` header.

    :return: JSON response with annotations.
    """

    return await cell_operations_service.annotate_adata_file(
        user=request_user,
        file=file.file,
        model_name=model_name,
        consensus_strategy=ConsensusStrategyType.ONTOLOGY_AWARE,
        normalize=normalize,
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
