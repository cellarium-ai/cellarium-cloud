import multiprocessing
import time
import typing as t

import sentry_sdk
import uvicorn
import uvicorn.config
from fastapi import APIRouter, FastAPI, Request, Response
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware, DispatchFunction
from starlette_context import context, plugins
from starlette_context.middleware import RawContextMiddleware
from starlette_context.plugins import Plugin

from casp.services import settings
from casp.services.api import exception_handlers
from casp.services.api.services import exceptions
from casp.services.logging import setup_logging


class RouterDef(t.NamedTuple):
    router: APIRouter
    prefix: str = "/api"
    tags: t.List[str] = []


class ExceptionHandlerDef(t.NamedTuple):
    exception: Exception
    handler: t.Callable


class CloudTraceContextPlugin(Plugin):
    key = "X-Cloud-Trace-Context"


BASE_PLUGGINS: t.Sequence[Plugin] = (
    # Adds an x-request-id header to all responses and makes it available to the request context.
    plugins.RequestIdPlugin(),
    # Adds an x-correlation-id header to all responses and makes it available to the request context.
    plugins.CorrelationIdPlugin(),
    # Extracts the user-agent header and makes it available to the request context.
    plugins.UserAgentPlugin(),
    # Extracts the x-cloud-trace-context header provided by CloudRun and makes it available to the request context.
)


BASE_ERROR_HANDLERS: t.List[ExceptionHandlerDef] = [
    ExceptionHandlerDef(exception=exceptions.AccessDeniedError, handler=exception_handlers.access_denied_error_handler),
    ExceptionHandlerDef(exception=exceptions.InvalidInputError, handler=exception_handlers.invalid_input_error_handler),
]


class TimingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next: DispatchFunction) -> Response:
        start_time = time.time()
        response = await call_next(request)
        request_time = time.time() - start_time
        context["request_time"] = request_time
        return response


class CASService(FastAPI):
    """
    Wrapper around FastAPI that provides some common configuration for CAS services.
    """

    def __init__(
        self,
        port: int,
        plugins: t.Optional[t.Sequence[Plugin]] = None,
        routers: t.List[RouterDef] = None,
        exception_handlers: t.Optional[t.List[ExceptionHandlerDef]] = None,
        sentry_application_id: str = "",
        *args,
        **kwargs,
    ):
        """
        Initialize the service with the given configuration.

        @param port: The port to use when running the service.
        @param plugins: A list of middleware plugins to use for the service.
        @param routers: A list of routers to include in the service.
        @param exception_handlers: A list of exception handlers to include in the service in addition to the default ones
                                   as defined in BASE_ERROR_HANDLERS
        @param sentry_application_id: The application ID to use when reporting errors to Sentry.
        """

        self.port = port

        # Configure Sentry integration
        sentry_sdk.init(
            dsn=settings.SENTRY_DSN,
            server_name=sentry_application_id,
            enable_tracing=settings.SENTRY_ENABLE_TRACING,
            profiles_sample_rate=settings.SENTRY_PROFILES_SAMPLE_RATE,
            traces_sample_rate=settings.SENTRY_TRACES_SAMPLE_RATE,
        )

        # # Init stackdriver logging
        # client = google.cloud.logging.Client()
        # client.setup_logging()

        setup_logging(settings.LOG_LEVEL, json_logs=settings.ENVIRONMENT != "local")

        # Configure middleware
        if plugins:
            _plugins = BASE_PLUGGINS + plugins
        else:
            _plugins = BASE_PLUGGINS
        middleware = [Middleware(RawContextMiddleware, plugins=_plugins), Middleware(TimingMiddleware)]

        # Perform basic initialization
        super().__init__(
            *args,
            **kwargs,
            version=settings.APP_VERSION,
            docs_url="/api/docs",
            redoc_url="/api/redoc",
            swagger_ui_parameters={"displayRequestDuration": True},
            middleware=middleware,
        )

        # Register routers
        for router_def in routers:
            self.include_router(router=router_def.router, prefix=router_def.prefix, tags=router_def.tags)

        # Register exception handlers
        for handler_def in exception_handlers or []:
            self.add_exception_handler(handler_def.exception, handler_def.handler)

    def run(self):
        """
        Launch the service using Uvicorn.
        """
        num_workers = 2 if settings.ENVIRONMENT == "local" else multiprocessing.cpu_count() * 2 + 1

        uvicorn.run(
            "main:application",
            host=settings.DEFAULT_SERVICE_HOST,
            port=self.port,
            workers=num_workers,
            log_level=settings.LOG_LEVEL,
            log_config=settings.LOG_CONFIG,
        )
