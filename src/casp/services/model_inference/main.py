from casp.services import settings
from casp.services.app import CASService, RouterDef
from casp.services.model_inference.routers import model_embed_router

application = CASService(
    title="Cellarium Cloud Model Inference",
    description="Cellarium Cloud Model Inference API",
    plugins=None,
    routers=[
        RouterDef(router=model_embed_router, tags=["model-inference"]),
    ],
    sentry_application_id="model-inference-service",
    local_port=settings.MODEL_SERVICE_PORT,
)

if __name__ == "__main__":
    application.run()
