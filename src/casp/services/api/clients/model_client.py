import asyncio
import json
import logging
import typing as t

import aiohttp
from aiohttp import client_exceptions
from starlette_context import context

from casp.services import settings
from casp.services.api.clients import exceptions
from casp.services.constants import HeaderKeys

JSON = t.Union[t.Dict[str, t.Any], t.List[t.Any]]

logger = logging.getLogger(__name__)


class HTTPAsyncClient:
    @staticmethod
    def raise_response_exception(status_code: int, detail: str) -> None:
        """
        Raise an exception based on the status code returned by the server, including the detail message

        :param status_code: HTTP status code
        :param detail: Detail message returned by the server
        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        message = f"Server returned status code {status_code}, Detail: {detail}"
        if status_code == 401:
            raise exceptions.HTTPError401(message)
        elif status_code == 403:
            raise exceptions.HTTPError403(message)
        elif status_code == 500 or status_code == 502 or status_code == 503 or status_code == 504:
            raise exceptions.HTTPError5XX(message)
        else:
            raise exceptions.HTTPResponseStatusCodeError(message)

    @classmethod
    async def _aiohttp_async_post(
        cls, url: str, form_data: t.Optional[aiohttp.FormData] = None, headers: t.Optional[t.Dict[str, t.Any]] = None
    ) -> JSON:
        """
        Make an async POST request to backend service with timeout and SSL context, handle exceptions and return JSON.

        Create custom SSL context with `certifi` package to avoid SSL errors. Use :class:`aiohttp.ClientTimeout` to set
        timeout for the request. Handle exceptions and raise custom exceptions.

        :param url: Endpoint URL
        :param form_data: :class:`aiohttp.FormData` object
        :param headers: Headers to include in the request

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError, HTTPClientError

        :return: JSON response
        """
        timeout = aiohttp.ClientTimeout(
            total=settings.AIOHTTP_CLIENT_TOTAL_TIMEOUT_SECONDS, sock_read=settings.AIOHTTP_CLIENT_TOTAL_TIMEOUT_SECONDS
        )

        headers = {} if headers is None else headers
        # Forward along the trace id if it exists
        if HeaderKeys.cloud_trace_context in context and context[HeaderKeys.cloud_trace_context] is not None:
            headers[HeaderKeys.cloud_trace_context.value] = context[HeaderKeys.cloud_trace_context]
        # Forward along the user auth if it exists
        if HeaderKeys.authorization in context and context[HeaderKeys.authorization] is not None:
            headers[HeaderKeys.authorization.value] = context[HeaderKeys.authorization]

        async with aiohttp.ClientSession(timeout=timeout) as session:
            try:
                async with session.post(url, data=form_data, headers=headers) as response:
                    status_code = response.status
                    if status_code < 200 or status_code >= 300:
                        try:
                            response_body = await response.json()
                            response_detail = response_body["detail"]
                        except (json.decoder.JSONDecodeError, client_exceptions.ClientResponseError):
                            response_detail = await response.text()
                        except KeyError:
                            logger.info("Response body doesn't have a 'detail' key, returning full response body")
                            response_detail = str(await response.json())

                        cls.raise_response_exception(status_code=status_code, detail=response_detail)
                    try:
                        return await response.json()
                    except client_exceptions.ClientPayloadError as e:
                        raise exceptions.ClientError(f"Failed to parse response body: {e.__class__.__name__}: {e}")

            except (client_exceptions.ServerTimeoutError, asyncio.TimeoutError) as e:
                raise exceptions.ClientTimeoutError(f"Client Timeout Error: {e}")
            except (client_exceptions.ClientConnectionError, client_exceptions.ClientResponseError) as e:
                raise exceptions.ClientError(f"Client Connection Error: {e.__class__.__name__}: {e}")
            except client_exceptions.ClientError as e:
                raise exceptions.ClientError(f"Unknown Error: {e.__class__.__name__}: {e}")

    @classmethod
    async def post(
        cls,
        url: str,
        data: t.Optional[t.Dict[str, t.Any]] = None,
        files: t.Optional[t.List[t.Dict[str, t.Any]]] = None,
        headers: t.Optional[t.Dict[str, t.Any]] = None,
    ) -> JSON:
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
        return await cls._aiohttp_async_post(url=url, form_data=form_data, headers=headers)


class ModelInferenceClient(HTTPAsyncClient):
    BACKEND_URL: str = settings.MODEL_SERVER_URL

    @classmethod
    def _get_endpoint_url(cls, endpoint: str) -> str:
        """
        Configure a specific method endpoint from backend url and endpoint

        :param endpoint: Endpoint string without a leading slash
        :return: Full url with backend domains/subdomains and endpoint joint
        """
        return f"{cls.BACKEND_URL}/{endpoint}"

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
        return await cls.post(url=url_endpoint, data=request_data, files=request_files)
