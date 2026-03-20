Changelog
=========

All notable changes to Cellarium Cloud will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.


1.8.0 - 2026-03-20
------------------

Added
~~~~~
- Added ``pyproject.toml`` and Poetry-based dependency management for the project

Changed
~~~~~~~
- Restructured the codebase from the legacy ``casp`` layout into ``cellarium.cas_backend`` with separate
  ``apps`` and ``core`` modules
- Renamed the main API service to CAS Compute and updated related configuration, documentation, and templates
- Upgraded the runtime and container stack to Python 3.12, simplified Docker base images, and refreshed entrypoints
- Reworked GitHub Actions and deployment workflows around service-specific build and deploy inputs
- Replaced the previous formatting and linting toolchain with Ruff and moved away from the older tox-driven flow

Fixed
~~~~~
- Fixed environment loading, import paths, and SQL template locations after the package restructuring
- Resolved remaining refactor regressions in tests, Cloud Run deployment wiring, and floating-point assertions


1.7.3 - 2025-03-20
------------------

Changed
~~~~~~~
- Updated the default Cloud Run deployment configuration
- Optimized ``Dockerfile.ns.pytorch`` to build a faster deployment image


1.7.2 - 2025-03-19
------------------

Added
~~~~~
- Added dedicated unit tests for the reimplemented summary statistics strategy

Changed
~~~~~~~
- Reimplemented the summary statistics consensus strategy used by cell operations
- Refactored API service construction to use dependency injection consistently
- Refreshed the unit-test lifecycle and SQLite-backed test database setup

Fixed
~~~~~
- Locked ``starlette`` and updated ``aiohttp`` and the PyTorch Docker image dependencies to restore compatibility


1.7.1 - 2024-11-13
------------------

Added
~~~~~
- Added more filtering criteria to SQL-template-driven BigQuery extraction workflows

Changed
~~~~~~~
- Refactored the ``bq_ops`` scripts and refreshed related documentation
- Updated dependencies to keep the audit workflow passing

Fixed
~~~~~
- Fixed non-random sampling behavior in extract workflows
- Fixed a timestamp handling bug affecting data processing


1.7.0 - 2024-10-25
------------------

Changed
~~~~~~~
- Lowered the initial user quota defaults while allowing quota increases through the REST API
- New users now default to non-admin accounts when created through the admin flow
- Updated the release version and minimum supported client version
- Refreshed the welcome email templates sent to new users


1.6.2 - 2024-10-21
------------------

Fixed
~~~~~
- Pinned the ``WTForms`` dependency to restore reliable admin installs


1.6.1 - 2024-10-21
------------------

Changed
~~~~~~~
- Updated the welcome and new-key email templates

Fixed
~~~~~
- User activity totals now default to ``0`` instead of ``None`` for processed cells and requests


1.6.0 - 2024-10-09
------------------

Added
~~~~~
- Added lifetime cell quotas, including database support and service logic for tracking them

Changed
~~~~~~~
- Raised the minimum supported client version for the lifetime-quota release


1.4.7 - 2024-09-06
------------------

Added
~~~~~
- Added a client compatibility endpoint so clients can check whether their version is supported
- Added a script for building the cell type ontology resources used by the consensus engine
- Added ``pip-audit`` checks to CI to catch vulnerable dependencies earlier

Changed
~~~~~~~
- Merged the model inference service into the API service and cached model assets to reduce repeated downloads
- Updated onboarding emails and raised the default weekly quota to 100K cells
- Allowed each deployment to set its own Sentry environment and refreshed release documentation

Fixed
~~~~~
- Reduced Sentry noise from infrastructure-level failures and pinned dependency versions with breaking changes
- Added stricter annotate input validation and removed redundant AnnData reads in the request flow


1.4.6 - 2024-07-03
------------------

Added
~~~~~
- Added user feedback capture, including redirect and opt-out support
- Added support for propagating ``x-client-session-id`` through request context and logs
- Added an integration-test workflow that boots a test server and exercises the API end to end
- Added description fields for models and indexes in the database schema

Changed
~~~~~~~
- Made deployed index IDs optional to support more flexible index configuration
- Updated the cell-related models and removed token display from the admin UI
- Updated deployment labels, email defaults, and Cloud SQL connector configuration

Fixed
~~~~~
- Resolved the multiple-head Alembic migration issue
- Fixed admin database connection management and user access logging


1.4.5 - 2024-06-04
------------------

Added
~~~~~
- Added opaque API keys with admin flows for creating, emailing, and authenticating user keys
- Added centralized logging and tracing support across services and downstream model calls
- Added weekly user cell quotas and quota-aware API behavior
- Added GitHub Actions support for building and deploying Cloud Run services
- Added support for connecting to PostgreSQL over VPC in cloud environments

Changed
~~~~~~~
- Centralized application configuration and generalized API exception handling
- Moved the development ``.env`` file under ``src/settings`` and updated deployment and secrets documentation
- Reworked Docker builds so secrets are no longer baked into images

Fixed
~~~~~
- Fixed logging credential loading and Sentry context setup order
- Fixed email template loading and related authentication edge cases


1.4.4 - 2024-05-02
------------------

Added
~~~~~
- Added REST API support for vector search indexes through a dedicated matching client layer
- Added request and cells-processed activity tracking in a dedicated user activity table
- Added the first consensus engine implementation, including ontology-aware strategies
- Added cell metadata and index-related database models plus supporting admin tooling

Changed
~~~~~~~
- Switched index search calls to an async REST client and stopped requesting unnecessary feature vectors
- Expanded tests and coverage around matching, consensus, and cell operations behavior

Fixed
~~~~~
- Fixed repeated matching-client initialization during search requests
- Fixed an ontology-aware strategy failure caused by missing temporary tables





1.4.2 - 2024-02-13
----------------------

Added
~~~~~
- Added Documentation in rst format for Sphinx build in `docs/` directory

Changed
~~~~~~~
- ``src/casp/services/README.md`` is now ``casp/services/README.rst``
- ``src/casp/services/admin/README.md`` is now ``casp/services/admin/README.rst``
- ``src/casp/services/api/README.md`` is now ``casp/services/api/README.rst``
- ``src/casp/services/model_inference/README.md`` is now ``casp/services/model_inference/README.rst``
- ``src/casp/services/db/README.md`` is now ``casp/services/db/README.rst``
- ``src/casp/services/wdl_workflows/README.md`` is now ``casp/services/wdl_workflows/README.rst``
- ``docs/VERTEX_AI_MATCHING_ENGINE.md`` is now ``docs/vertex_ai_matching_engine.rst``
- Requests to matching engine are now made in batches of 5 to avoid overloading Vertex API
- Requests to matching engine are now made with retry logic to account for errors caused by rapid increases in traffic
