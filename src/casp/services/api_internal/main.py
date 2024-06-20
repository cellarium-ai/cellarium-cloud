from casp.services import settings
from casp.services.api_internal.routers import api_internal_router
from casp.services.app import CASService, RouterDef

application = CASService(
    title="Cellarium Cloud Internal API",
    description="Cellarium Cloud Application API for internal use",
    plugins=None,
    routers=[
        RouterDef(router=api_internal_router, tags=["ml-service-management"]),
    ],
    sentry_application_id="api-internal-service",
    port=settings.API_INTERNAL_SERVICE_PORT,
)

if __name__ == "__main__":
    application.run()
