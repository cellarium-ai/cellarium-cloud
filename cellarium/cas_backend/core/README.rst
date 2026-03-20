CAS Backend Core
================

The ``cellarium.cas_backend.core`` package contains shared runtime infrastructure used by the deployable apps.
It provides configuration loading, FastAPI service wiring, authentication helpers, database access, data managers,
logging utilities, and other common building blocks.

What Lives Here
---------------

- ``app.py``: shared ``CASService`` wrapper and middleware setup for FastAPI services
- ``config.py``: environment-driven settings and repository path resolution
- ``auth/``: JWT and opaque token helpers shared across services
- ``data_managers/``: BigQuery-oriented query composition and data access code
- ``db/``: SQLAlchemy models, migration config, and persistence helpers
- ``utils/``: shared integration utilities such as email and Google Cloud helpers

Public Entry Points
-------------------

- ``cellarium.cas_backend.core.app.CASService`` is the common application wrapper for service processes.
- ``cellarium.cas_backend.core.config.settings`` is the shared settings object used throughout the codebase.
- ``cellarium.cas_backend.core.db`` exposes the relational data layer used by admin and compute services.
- ``cellarium.cas_backend.core.data_managers`` exposes the warehouse-facing data access layer used by compute flows.

Dependencies
------------

The core package is consumed primarily by:

- ``cellarium.cas_backend.apps.compute``
- ``cellarium.cas_backend.apps.admin``
- internal scripts and deployment entrypoints

Core depends on shared third-party infrastructure such as PostgreSQL, BigQuery, Cloud Run configuration, and optional
Sentry integration depending on the environment.

Local Development Notes
-----------------------

Configuration is loaded from ``settings/.env`` at the repository root. Copy ``settings/sample_env`` to
``settings/.env`` before running services locally.

Related Documentation
---------------------

- ``cellarium/cas_backend/core/db/README.rst``
- ``cellarium/cas_backend/core/data_managers/README.rst``
- ``docs/source/modules/secrets.rst``
