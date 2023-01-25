import typing as t

import aiohttp

from casp.services import settings


class CASAPIAsyncClient:
    @classmethod
    async def post(cls, url: str, data: t.Dict):
        async with aiohttp.ClientSession() as session:
            async with session.post(url, data=data) as resp:
                return await resp.text()

    @classmethod
    async def call_model_service(cls, file_to_embed):
        return await cls.post(url=settings.MODEL_SERVER_URL, data={"file": file_to_embed})
