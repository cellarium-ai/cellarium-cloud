from cellarium.cas_backend.apps.compute.exceptions import APIInternalError


class VectorSearchConfigurationError(APIInternalError):
    http_code: int = APIInternalError.http_code
