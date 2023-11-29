import multiprocessing

import uvicorn
from fastapi import FastAPI

from casp.services.model_inference import routers

application = FastAPI()
application.include_router(router=routers.model_embed_router, prefix="/api/model_inference", tags=["Model Inference"])

if __name__ == "__main__":
    # Run model server
    uvicorn.run("main:application", host="0.0.0.0", port=8000, workers=multiprocessing.cpu_count())
