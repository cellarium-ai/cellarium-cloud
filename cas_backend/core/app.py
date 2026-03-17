import logging
import multiprocessing
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
from casp.services._auth.exceptions import TokenException
from casp.services.api import exception_handlers
from casp.services.api.data_manager.exceptions import NotFound
from casp.services.api.services.exceptions import APIBaseException
from casp.services.constants import ContextKeys, HeaderKeys, SentryTags

logger = logging.getLogger(__name__)


class RouterDef(t.NamedTuple):
    """
    Used to define a router to include in the service.
    """

    router: APIRouter
    prefix: str = "/api"
    tags: t.List[str] = []


class ExceptionHandlerDef(t.NamedTuple):
    """
    Used to define an exception handler for a specific exception type.
    """

    exception: Exception
    handler: t.Callable


class CloudTraceContextPlugin(Plugin):
    key = HeaderKeys.cloud_trace_context


class AuthorizationPlugin(Plugin):
    key = HeaderKeys.authorization


class ClientSessionPlugin(Plugin):
    key = HeaderKeys.client_session_id


class ClientActionPlugin(Plugin):
    key = HeaderKeys.client_action_id


BASE_PLUGGINS: t.Sequence[Plugin] = (
    # Extracts the user-agent header and makes it available to the request context.
    plugins.UserAgentPlugin(),
    # Extracts the x-cloud-trace-context header provided by CloudRun and makes it available to the request context.
    CloudTraceContextPlugin(),
    # Extracts the authorization header and makes it available to the request context.
    AuthorizationPlugin(),
    # Extracts the x-client-session-id header and makes it available to the request context.
    ClientSessionPlugin(),
    # Extracts the x-client-action-id header and makes it available to the request context.
    ClientActionPlugin(),
)


BASE_ERROR_HANDLERS: t.List[ExceptionHandlerDef] = [
    # Default error handler for APIBaseException.  Handles most API-thrown exceptions.
    ExceptionHandlerDef(exception=APIBaseException, handler=exception_handlers.api_base_exception_handler),
    # Special error handler for data_mananger not found exceptions in case they are not caught.
    ExceptionHandlerDef(exception=NotFound, handler=exception_handlers.not_found_error_handler),
    # Error handler for token-related exceptions.
    ExceptionHandlerDef(exception=TokenException, handler=exception_handlers.token_exception_handler),
    # Global error handler for all other exceptions.
    ExceptionHandlerDef(exception=Exception, handler=exception_handlers.global_error_handler),
]


class RequestClientMiddleware(BaseHTTPMiddleware):
    """
    Adds the request's client object to the request context.  This doesn't come from a
    header which is why this can't be a plugin.
    """

    async def dispatch(self, request: Request, call_next: DispatchFunction) -> Response:
        context[ContextKeys.client] = request.client
        response = await call_next(request)
        return response


class ContextToSentryMiddleware(BaseHTTPMiddleware):
    """
    Adds the value from the request context to the Sentry context.
    """

    def __init__(self, context_key: HeaderKeys, sentry_tag: SentryTags, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.context_key = context_key
        self.sentry_tag_name = sentry_tag

    async def dispatch(self, request: Request, call_next: DispatchFunction) -> Response:
        value = None
        if self.context_key in context:
            value = context[self.context_key]
            sentry_sdk.set_tag(self.sentry_tag_name, value)
        response = await call_next(request)
        if hasattr(response, "headers") and value is not None:
            response.headers[self.context_key] = value
        return response


class ClientSessionIdMiddleware(ContextToSentryMiddleware):
    """
    Adds the request's session id to the Sentry context.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(
            context_key=HeaderKeys.client_session_id,
            sentry_tag=SentryTags.client_session_id,
            *args,
            **kwargs,
        )


class ClientActionIdMiddleware(ContextToSentryMiddleware):
    """
    Adds the request's action id to the Sentry context.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(
            context_key=HeaderKeys.client_action_id,
            sentry_tag=SentryTags.client_action_id,
            *args,
            **kwargs,
        )


class SentryResponseDecoratorMiddleware(BaseHTTPMiddleware):
    """
    Adds the sentry traceid to the request response headers for support purposes.
    """

    async def dispatch(self, request: Request, call_next: DispatchFunction) -> Response:
        # Format is [trace_id]-[span_id]-[sampled].  We're only really interested in the trace_id.
        trace_id = sentry_sdk.get_traceparent()
        if trace_id:
            trace_id = trace_id.split("-")[0]
        context[ContextKeys.sentry_trace_id] = trace_id

        # Process the request
        response = await call_next(request)

        # This should always be true but adding check to be safe
        if hasattr(response, "headers"):
            response.headers[HeaderKeys.trace_id] = trace_id
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
        sentry_application_id: t.Optional[str] = None,
        *args,
        **kwargs,
    ):
        """
        Initialize the service with the given configuration.

        :param port: The port to use when running the service.
        :param plugins: A list of middleware plugins to use for the service.
        :param routers: A list of routers to include in the service.
        :param exception_handlers: A list of exception handlers to include in the service in addition to the default ones
                                   as defined in BASE_ERROR_HANDLERS
        :param sentry_application_id: The application ID to use when reporting errors to Sentry. If left as None, Sentry
                                      integration will not be enabled.
        """

        self.port = port

        # Configure Sentry integration
        if (
            sentry_application_id is not None
            and settings.SENTRY_DSN is not None
            and settings.SENTRY_ENVIRONMENT is not None
        ):
            sentry_sdk.init(
                dsn=settings.SENTRY_DSN,
                environment=settings.SENTRY_ENVIRONMENT,
                server_name=sentry_application_id,
                enable_tracing=settings.SENTRY_ENABLE_TRACING,
                traces_sample_rate=settings.SENTRY_TRACES_SAMPLE_RATE,
                profiles_sample_rate=settings.SENTRY_PROFILES_SAMPLE_RATE,
            )

        # Configure middleware
        if plugins:
            _plugins = BASE_PLUGGINS + plugins
        else:
            _plugins = BASE_PLUGGINS
        middleware = [
            Middleware(RawContextMiddleware, plugins=_plugins),
            Middleware(RequestClientMiddleware),
            Middleware(SentryResponseDecoratorMiddleware),
            Middleware(ClientSessionIdMiddleware),
            Middleware(ClientActionIdMiddleware),
        ]

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
        for handler_def in (exception_handlers or []) + BASE_ERROR_HANDLERS:
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
