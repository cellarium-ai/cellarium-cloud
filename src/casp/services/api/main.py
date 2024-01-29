import multiprocessing

import sentry_sdk
import uvicorn
from fastapi import FastAPI

from casp.services import settings
from casp.services.api import exception_handlers
from casp.services.api.routers import cell_operations_router, cellarium_general_router
from casp.services.api.services import exceptions

sentry_sdk.init(
    dsn=settings.SENTRY_DSN,
    enable_tracing=settings.SENTRY_ENABLE_TRACING,
    profiles_sample_rate=settings.SENTRY_PROFILES_SAMPLE_RATE,
    traces_sample_rate=settings.SENTRY_TRACES_SAMPLE_RATE,
)

application = FastAPI(
    title="Cellarium Cloud API",
    description="Cellarium Cloud Application API",
    version=settings.APP_VERSION,
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)
API_SERVICE_PREFIX = "/api"

application.include_router(router=cellarium_general_router, prefix=API_SERVICE_PREFIX, tags=["cellarium-general"])
application.include_router(router=cell_operations_router, prefix=API_SERVICE_PREFIX, tags=["cell-operations"])

application.add_exception_handler(exceptions.AccessDeniedError, exception_handlers.access_denied_error_handler)
application.add_exception_handler(exceptions.InvalidInputError, exception_handlers.invalid_input_error_handler)

if __name__ == "__main__":
    num_workers = multiprocessing.cpu_count() * 2 + 1
    uvicorn.run("main:application", host=settings.SERVER_HOST, port=settings.SERVER_PORT, workers=num_workers)
