import asyncio
import json
import typing as t

import aiohttp
import requests
from aiohttp import client_exceptions

from casp.clients import exceptions
from casp.services import settings

JSON = t.Union[t.Dict[str, t.Any], t.List[t.Any]]


class HTTPClient:
    BACKEND_URL: str

    @classmethod
    def _get_endpoint_url(cls, endpoint: str) -> str:
        """
        Configure a specific method endpoint from backend url and endpoint

        :param endpoint: Endpoint string without a leading slash
        :return: Full url with backend domains/subdomains and endpoint joint
        """
        return f"{cls.BACKEND_URL}/{endpoint}"

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

    def __validate_requests_response(self, response: requests.Response) -> None:
        """
        Validate requests response and raise an exception if response status code is not 200

        :param response: Response object

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        status_code = response.status_code
        if status_code < 200 or status_code >= 300:
            try:
                response_detail = response.json()["detail"]
            except (json.decoder.JSONDecodeError, requests.exceptions.JSONDecodeError, KeyError):
                response_detail = response.text

            self.raise_response_exception(status_code=status_code, detail=response_detail)

    def get(self, endpoint: str) -> requests.Response:
        """
        Make a GET request to backend service

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: Response object
        """
        url = self._get_endpoint_url(endpoint)
        response = requests.get(url=url)
        self.__validate_requests_response(response=response)
        return response

    def post(self, endpoint: str, data: t.Optional[t.Dict] = None) -> requests.Response:
        """
        Make a synchronous POST request to backend service.

        :param endpoint: Endpoint string without a leading slash
        :param data: Dictionary to include to HTTP POST request body

        :return: Response object
        """
        url = self._get_endpoint_url(endpoint)
        response = requests.post(url=url, json=data)
        self.__validate_requests_response(response=response)
        return response

    def get_json(self, endpoint: str) -> JSON:
        """
        Make a GET request to backend service and return JSON response

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        return self.get(endpoint=endpoint).json()

    def post_json(self, endpoint: str, data: t.Optional[t.Dict] = None) -> JSON:
        """
        Make a POST request to backend service and return JSON response

        :param endpoint: Endpoint string without a leading slash
        :param data: Dictionary to include to HTTP POST request body

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        return self.post(endpoint=endpoint, data=data).json()

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
                            print("Response body doesn't have a 'detail' key, returning full response body")
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
    async def async_post(
        cls, url: str, data: t.Optional[t.Dict[str, t.Any]] = None, files: t.Optional[t.List[t.Dict[str, t.Any]]] = None
    ) -> JSON:
        """
        Make an async POST request to backend service with timeout and SSL context, handle exceptions and return JSON.

        :param url: Endpoint URL
        :param data: Request body data
        :param files: List of dictionaries with file objects, filename and form data field name

        :return: JSON response
        """
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

        return await cls._aiohttp_async_post(url=url, form_data=form_data)
