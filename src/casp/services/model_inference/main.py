import multiprocessing

import sentry_sdk
import uvicorn
from fastapi import FastAPI

from casp.services import settings
from casp.services.model_inference import routers

sentry_sdk.init(
    dsn=settings.SENTRY_DSN,
    enable_tracing=settings.SENTRY_ENABLE_TRACING,
    profiles_sample_rate=settings.SENTRY_PROFILES_SAMPLE_RATE,
    traces_sample_rate=settings.SENTRY_TRACES_SAMPLE_RATE,
)

application = FastAPI(
    title="Cellarium Cloud Model Inference",
    description="Cellarium Cloud Model Inference API Documentation",
    version=settings.APP_VERSION,
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)
application.include_router(router=routers.model_embed_router, prefix="/api", tags=["model-inference"])

if __name__ == "__main__":
    num_workers = 2 if settings.ENVIRONMENT == "local" else multiprocessing.cpu_count()
    port = settings.MODEL_SERVICE_PORT if settings.ENVIRONMENT == "local" else settings.DEFAULT_SERVICE_PORT
    uvicorn.run("main:application", host=settings.DEFAULT_SERVICE_HOST, port=port, workers=num_workers)
