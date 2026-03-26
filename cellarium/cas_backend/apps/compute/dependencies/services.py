from cellarium.cas_backend.apps.compute.services.cell_operations_service import CellOperationsService
from cellarium.cas_backend.apps.compute.services.cell_quota_service import CellQuotaService
from cellarium.cas_backend.apps.compute.services.cellarium_general_service import CellariumGeneralService


def get_cell_operations_service() -> CellOperationsService:
    return CellOperationsService()


def get_cellarium_general_service() -> CellariumGeneralService:
    return CellariumGeneralService()


def get_cell_quota_service() -> CellQuotaService:
    return CellQuotaService()
