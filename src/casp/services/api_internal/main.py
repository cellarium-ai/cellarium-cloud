import multiprocessing

import uvicorn
from fastapi import FastAPI

from casp.services import settings
from casp.services.api_internal import exception_handlers, exceptions, routers

application = FastAPI(
    title="Cellarium Cloud API Internal",
    description="Cellarium Cloud Application API for internal services use only.",
    version=settings.APP_VERSION,
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)
API_SERVICE_PREFIX = "/api"

application.include_router(
    router=routers.api_internal_router, prefix=API_SERVICE_PREFIX, tags=["cellarium-api-internal"]
)

application.add_exception_handler(exceptions.DatabaseError, exception_handlers.database_unique_constraint_error_handler)

if __name__ == "__main__":
    num_workers = multiprocessing.cpu_count() * 2 + 1
    uvicorn.run("main:application", host=settings.SERVER_HOST, port=settings.SERVER_PORT, workers=num_workers)
