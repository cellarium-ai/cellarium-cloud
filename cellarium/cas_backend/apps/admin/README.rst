Admin Service
=============

The admin service is the internal Flask application used to manage relational data and operational state for CAS
Backend. It exposes administrative views for users, keys, models, indexes, and related database-backed records.

What Lives Here
---------------

- ``server.py``: local execution entrypoint
- ``__init__.py``: Flask app creation and extension setup
- ``views.py``: Flask-Admin views and admin workflows
- ``templates/``: admin UI and email templates
- ``entrypoint.sh``: container startup command

Public Entrypoints
------------------

- Local process: ``python -m cellarium.cas_backend.apps.admin.server``
- Container entrypoint: ``cellarium/cas_backend/apps/admin/entrypoint.sh``

Dependencies
------------

The admin service depends on:

- ``cellarium.cas_backend.core.config`` for environment settings
- ``cellarium.cas_backend.core.db`` for models, sessions, and migrations
- ``cellarium.cas_backend.core.utils.email_utils`` for admin-triggered emails
- ``settings/.env`` for runtime secrets

Local Development
-----------------

Install dependencies and ensure ``settings/.env`` is present, then run:

.. code-block:: shell

    poetry run python -m cellarium.cas_backend.apps.admin.server

By default the app serves the admin UI locally on Flask's default port. Authentication and database connectivity still
depend on the configured environment variables.

Deployment Notes
----------------

The admin image is built from ``deploy/docker/Dockerfile.admin`` and deployed through the shared Cloud Run workflows
and flavor configuration under ``deploy/cloudrun/``.

Related Documentation
---------------------

- ``cellarium/cas_backend/core/db/README.rst``
- ``docs/source/modules/deploying_in_cloud.rst``
- ``docs/source/modules/secrets.rst``
