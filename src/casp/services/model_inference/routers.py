from fastapi import APIRouter, File, Form, UploadFile

from casp.services.model_inference import schemas, services

model_embed_router = APIRouter()

model_inference_service = services.ModelInferenceService()


@model_embed_router.post("/embed", response_model=schemas.ModelEmbeddings)
def embed(file: UploadFile = File(), model_name: str = Form()) -> schemas.ModelEmbeddings:
    """
    Embed adata file using a specific model using Cellarium-ML model and pytorch

    :param file: File object of :class:`anndata.AnnData` object to embed.
    :param model_name: Model name to use for embedding.

    :return: ModelEmbeddings object
    """
    return model_inference_service.embed_adata_file(file_to_embed=file.file, model_name=model_name)
