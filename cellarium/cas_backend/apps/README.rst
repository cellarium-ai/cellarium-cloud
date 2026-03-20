CAS Backend Apps
================

The ``cellarium.cas_backend.apps`` package contains deployable services and other app-scoped code.
These packages own HTTP entrypoints, dependency wiring, request/response schemas, and domain services that sit on top
of the shared ``core`` package.

What Lives Here
---------------

- ``admin/``: Flask-based internal administration service
- ``compute/``: FastAPI-based compute-facing service and its domain logic
- ``model_inference/``: model-inference support code used by inference-oriented app flows

Public Entrypoints
------------------

- ``cellarium.cas_backend.apps.compute.main:application``
- ``cellarium.cas_backend.apps.admin.server``
- shell entrypoints under each deployable app for container startup

Dependencies
------------

App packages depend on ``cellarium.cas_backend.core`` for configuration, auth, shared service setup, data access, and
database integration. Operational deployment details live outside this package under ``deploy/``.

Documentation Boundaries
------------------------

Keep app READMEs focused on module-specific concerns:

- purpose and responsibilities
- important subpackages
- service entrypoints
- local run/test commands that are specific to that app

Cross-cutting runbooks such as CI/CD, secrets, Cloud Run deployment, and shared local setup belong in the dedicated
Sphinx pages under ``docs/source/modules/``.
