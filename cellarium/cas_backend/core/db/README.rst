Database Module
===============

The ``cellarium.cas_backend.core.db`` package contains the relational persistence layer for CAS Backend.
It defines SQLAlchemy models, Alembic migration configuration, and helper operations used by the app layer.

What Lives Here
---------------

- ``models/``: SQLAlchemy models for users, cell-related records, and ML management metadata
- ``ops.py``: shared database operations
- ``alembic.ini`` and ``migrations/``: Alembic configuration and migration history
- ``field_type_decorators.py``: shared SQLAlchemy field helpers

Public Entry Points
-------------------

- ``cellarium.cas_backend.core.db.models``
- ``cellarium.cas_backend.core.db.ops``
- Alembic config at ``cellarium/cas_backend/core/db/alembic.ini``

Dependencies
------------

This package is used by the admin and compute services. Runtime configuration comes from
``cellarium.cas_backend.core.config.settings``.

Environment and Connection
--------------------------

Local development reads database settings from ``settings/.env``. The code supports:

- local PostgreSQL connections via ``DB_HOST``, ``DB_PORT``, ``DB_USER``, ``DB_PASSWORD``, and ``DB_NAME``
- Cloud SQL connections via ``DB_INSTANCE_UNIX_SOCKET`` or ``DB_PRIVATE_IP`` with the same credential variables

The SQLAlchemy URI is assembled in ``cellarium.cas_backend.core.config``.

Database Migrations
-------------------

Generate a migration:

.. code-block:: shell

    alembic -c cellarium/cas_backend/core/db/alembic.ini revision --autogenerate -m "describe_change"

Apply migrations:

.. code-block:: shell

    alembic -c cellarium/cas_backend/core/db/alembic.ini upgrade head

Related Documentation
---------------------

- ``docs/source/modules/secrets.rst``
- ``docs/source/modules/running_locally.rst``
