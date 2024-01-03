import typing as t

from casp.clients import HTTPClient
from casp.services import settings

JSON = t.Union[t.Dict[str, t.Any], t.List[t.Any]]


class ModelInferenceClient(HTTPClient):
    BACKEND_URL: str = settings.MODEL_SERVER_URL

    @classmethod
    async def call_model_embed(cls, file_to_embed: bytes, model_name: str) -> JSON:
        """
        Call model microservice to embed adata file using a specific model from Cellarium Cloud infrastructure.

        :param file_to_embed: Instance :class:`anndata.AnnData` file to embed
        :param model_name: Model name to use for embedding

        :return: JSON response with base64 encoded embeddings and obs_ids
        """
        request_data = {"model_name": model_name}
        request_files = [{"file": file_to_embed, "filename": "adata.h5ad", "form_data_field_name": "file"}]
        url_endpoint = cls._get_endpoint_url(endpoint="api/embed")
        return await cls.async_post(url=url_endpoint, data=request_data, files=request_files)
