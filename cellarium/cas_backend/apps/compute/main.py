from cellarium.cas_backend.apps.compute.routers import cell_operations_router, cellarium_general_router
from cellarium.cas_backend.core.app import CASService, RouterDef
from cellarium.cas_backend.core.config import settings

application = CASService(
    title="CAS Backend API",
    description="CAS Backend Application API",
    plugins=None,
    routers=[
        RouterDef(router=cellarium_general_router, tags=["cellarium-general"]),
        RouterDef(router=cell_operations_router, tags=["cell-operations"]),
    ],
    sentry_application_id="api-service",
    port=settings.API_SERVICE_PORT,
    app_module_path="cellarium.cas_backend.apps.compute.main:application",
)

if __name__ == "__main__":
    application.run()
