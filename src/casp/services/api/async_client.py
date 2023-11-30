import typing as t

import aiohttp

from casp.services import settings
from casp.services.api import exceptions


AIOHTTP_TOTAL_TIMEOUT_SECONDS = 1200
AIOHTTP_READ_TIMEOUT_SECONDS = 1100


class CASAPIAsyncClient:
    @classmethod
    async def post(
        cls, url: str, data: t.Optional[t.Dict[str, t.Any]] = None, files: t.Optional[t.List[t.Dict[str, t.Any]]] = None
    ) -> str:
        form_data = aiohttp.FormData()

        if files is not None:
            for file_obj in files:
                file = file_obj["file"]
                filename = file_obj["filename"]
                form_data_field_name = file_obj["form_data_field_name"]
                form_data.add_field(form_data_field_name, file, filename=filename)

        if data is not None:
            for key, value in data.items():
                form_data.add_field(key, value)

        # Client Session Arguments
        timeout = aiohttp.ClientTimeout(total=AIOHTTP_TOTAL_TIMEOUT_SECONDS, sock_read=AIOHTTP_READ_TIMEOUT_SECONDS)

        async with aiohttp.ClientSession(timeout=timeout) as session:
            async with session.post(url, data=form_data) as response:
                status = response.status

                if status < 200 or status >= 300:
                    response_body = await response.text()
                    raise exceptions.ServiceAPIException(
                        f"Model Service returned status code {status}, Detail: {response_body}"
                    )
                return await response.text()

    @classmethod
    async def call_model_service(cls, file_to_embed: bytes, model_name: str) -> str:
        request_data = {"model_name": model_name}
        request_files = [{"file": file_to_embed, "filename": "adata.h5ad", "form_data_field_name": "file"}]
        return await cls.post(url=settings.MODEL_SERVER_URL, data=request_data, files=request_files)
