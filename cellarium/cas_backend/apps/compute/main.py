from cellarium.cas_backend.apps.compute.routers import cell_operations_router, cellarium_general_router
from cellarium.cas_backend.core.app import CASService, RouterDef
from cellarium.cas_backend.core.config import settings

application = CASService(
    title="Cellarium Cloud API",
    description="Cellarium Cloud Application API",
    plugins=None,
    routers=[
        RouterDef(router=cellarium_general_router, tags=["cellarium-general"]),
        RouterDef(router=cell_operations_router, tags=["cell-operations"]),
    ],
    sentry_application_id="api-service",
    port=settings.API_SERVICE_PORT,
)

if __name__ == "__main__":
    application.run()
