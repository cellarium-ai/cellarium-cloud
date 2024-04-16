from casp.services import settings
from casp.services.api.routers import cell_operations_router, cellarium_general_router
from casp.services.app import CASService, RouterDef

application = CASService(
    title="Cellarium Cloud API",
    description="Cellarium Cloud Application API",
    plugins=None,
    routers=[
        RouterDef(router=cellarium_general_router, tags=["cellarium-general"]),
        RouterDef(router=cell_operations_router, tags=["cell-operations"]),
    ],
    sentry_application_id="api-service",
    local_port=settings.API_SERVICE_PORT,
)

if __name__ == "__main__":
    application.run()
