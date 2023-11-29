import multiprocessing

import uvicorn
from fastapi import FastAPI

from casp.services import settings
from casp.services.api.routers import cell_analysis_router, cellarium_general_router

application = FastAPI(
    title="Cellarium Cloud API",
    description="Cellarium Cloud Application API",
    version=settings.APP_VERSION,
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)
API_SERVICE_PREFIX = "/api"

application.include_router(router=cellarium_general_router, prefix=API_SERVICE_PREFIX, tags=["Cellarium General"])
application.include_router(router=cell_analysis_router, prefix=API_SERVICE_PREFIX, tags=["Cell Analysis"])

from casp.services.api import exception_handlers  # noqa: E402, F401

if __name__ == "__main__":
    num_workers = multiprocessing.cpu_count() * 2 + 1
    uvicorn.run("main:application", host=settings.SERVER_HOST, port=settings.SERVER_PORT, workers=num_workers)
