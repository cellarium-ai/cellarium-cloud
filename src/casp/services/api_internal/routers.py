from fastapi import APIRouter

from casp.services.api_internal import schemas, services

api_internal_router = APIRouter()

service = services.EmbeddingModelRegistryService()


@api_internal_router.post(path="/cas-model", response_model=schemas.CASModelOut, status_code=201)
async def create_embedding_model(model: schemas.CASModelInCreate):
    return service.create_embedding_model(
        model_name=model.model_name,
        model_file_path=model.model_file_path,
        embedding_dimension=model.embedding_dimension,
        bq_dataset_name=model.bq_dataset_name,
        schema_name=model.schema_name,
    )


@api_internal_router.post(path="/cas-index", response_model=schemas.CASIndexOut, status_code=201)
async def create_index(index: schemas.CASIndexInCreate):
    return service.create_index(
        model_name=index.model_name,
        index_name=index.index_name,
        num_neighbors=index.num_neighbors,
        deployed_index_id=index.deployed_index_id,
        endpoint_id=index.endpoint_id,
        embedding_dimension=index.embedding_dimension,
    )


@api_internal_router.patch(path="cas-index/{index_name}", response_model=schemas.CASIndexOut, status_code=200)
async def update_index(index_name, item: schemas.CASIndexInUpdate):
    return service.update_index(index_name=index_name, index_schema_item=item)
