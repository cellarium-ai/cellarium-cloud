Authentication Helpers
======================

The ``cellarium.cas_backend.core.auth`` package contains shared authentication utilities used by deployable services.
It is not a standalone service and is intended to be imported by app-layer authentication dependencies.

What Lives Here
---------------

- ``jwt_token.py``: JWT helpers
- ``opaque_token.py``: database-backed opaque token helpers
- ``exceptions.py``: auth-related exceptions

Dependencies
------------

The auth helpers depend on the core database layer and shared configuration.
