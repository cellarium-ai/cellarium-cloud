import multiprocessing

import sentry_sdk
import uvicorn
from fastapi import FastAPI

from casp.services import settings
from casp.services.model_inference import routers

sentry_sdk.init(
    dsn=settings.SENTRY_DSN,
    server_name="model-inference-service",
    enable_tracing=settings.SENTRY_ENABLE_TRACING,
    profiles_sample_rate=settings.SENTRY_PROFILES_SAMPLE_RATE,
    traces_sample_rate=settings.SENTRY_TRACES_SAMPLE_RATE,
)

application = FastAPI()
application.include_router(router=routers.model_embed_router, prefix="/api/model_inference", tags=["model-inference"])

if __name__ == "__main__":
    # Run model server
    uvicorn.run("main:application", host="0.0.0.0", port=8000, workers=multiprocessing.cpu_count())
