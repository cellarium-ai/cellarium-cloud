Compute Service
===============

The compute service is the primary FastAPI application in CAS Backend. It serves authenticated API traffic, coordinates
annotation and query workflows, and composes shared core services, data managers, and external integrations such as
BigQuery, Cloud SQL, and Vertex AI Matching Engine.

What Lives Here
---------------

- ``main.py``: FastAPI application entrypoint
- ``routers/``: HTTP route definitions
- ``schemas/``: request and response models
- ``dependencies/``: FastAPI dependency providers, including auth and service wiring
- ``services/``: domain logic for compute workflows
- ``clients/``: clients for external compute-facing integrations
- ``vector_search/``: provider-specific vector search adapters used by compute
- ``entrypoint.sh``: container startup command

Public Entrypoints
------------------

- Local process: ``python -m cellarium.cas_backend.apps.compute.main``
- ASGI app: ``cellarium.cas_backend.apps.compute.main:application``
- Container entrypoint: ``cellarium/cas_backend/apps/compute/entrypoint.sh``

Dependencies
------------

The compute service depends on:

- ``cellarium.cas_backend.core.app`` for common FastAPI service wiring
- ``cellarium.cas_backend.core.auth`` for authentication helpers
- ``cellarium.cas_backend.core.data_managers`` for BigQuery and warehouse access
- ``cellarium.cas_backend.core.db`` for quota, user, and metadata persistence
- vector search backends such as Vertex AI Matching Engine and TileDB for production workloads

Local Development
-----------------

Create ``settings/.env`` and install dependencies, then run:

.. code-block:: shell

    poetry run python -m cellarium.cas_backend.apps.compute.main

Useful commands during development:

.. code-block:: shell

    make test
    poetry run pytest tests/unit -k compute

The service publishes OpenAPI docs at ``http://localhost:8000/api/docs`` when running locally.

Deployment Notes
----------------

The compute image is built from ``deploy/docker/Dockerfile.compute`` and deployed through the shared Cloud Run
workflows. Deployment flavor and Cloud Run sizing live under ``deploy/cloudrun/``.

Related Documentation
---------------------

- ``docs/source/modules/deploying_in_cloud.rst``
- ``docs/source/modules/secrets.rst``
