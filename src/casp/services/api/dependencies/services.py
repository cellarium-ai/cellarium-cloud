from casp.services.api import services


def get_cell_operations_service() -> services.CellOperationsService:
    return services.CellOperationsService()


def get_cellarium_general_service() -> services.CellariumGeneralService:
    return services.CellariumGeneralService()


def get_cell_quota_service() -> services.CellQuotaService:
    return services.CellQuotaService()
