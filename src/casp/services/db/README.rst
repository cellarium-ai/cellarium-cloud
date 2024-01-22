Database Service
================

This module represents the codebase for connecting to the SQL database and performing operations on it. It is not a deployable service,
but rather a library that is used by other services.

.. note:: This service requires a running database cluster.

Setting up a Database Cluster
-----------------------------

Create a db cluster:

.. code-block:: bash

    gcloud sql instances create cas-db-cluster \
    --database-version=POSTGRES_14 \
    --cpu=1 \
    --memory 3.75GB \
    --storage-size 10GB \
    --require-ssl \
    --region=us-central

Create a cas-db-user:

.. code-block:: bash

    SECRET_PASSWORD='secretPassW0rd'

    gcloud sql users create cas-db-user \
    --instance=cas-db-cluster \
    --password=$SECRET_PASSWORD

Create a cas-db:

.. code-block:: bash

    gcloud sql databases create cas-db \
    --instance=cas-db-cluster

Code Base Info
--------------

Database module ``db`` consist of:

- ``migrations`` Database migration history;
- ``models.py`` All database models;
- ``ops.py`` Data operations;
- ``alembic.ini`` Provides metadata for migration manager;

Database Adapter and ORM
------------------------

`SQLAlchemy <https://www.sqlalchemy.org/>`_ is used for db object-relational mapping.
`pg8000 <https://pypi.org/project/pg8000/>`_ is used as a database adapter.

Environment Variables and Connection
-------------------------------------

Locally, the connection is approached through a regular PostgreSQL connection (by port). In the cloud, it's done through a proxy. To connect to the proxy it uses a Unix socket.

``SQLALCHEMY_DATABASE_URI`` is used by `SQLAlchemy <https://www.sqlalchemy.org/>`_ to connect and configured from the following secret environment variables:

- ``DB_HOST``, ``DB_PORT``, ``DB_USER``, ``DB_PASSWORD``, ``DB_NAME`` in the local environment.
- ``DB_USER``, ``DB_PASSWORD``, ``DB_NAME``, ``DB_INSTANCE_UNIX_SOCKET`` in development and production environments.

Database Migrations
-------------------

`Alembic <https://alembic.sqlalchemy.org/en/latest/>`_ is used for managing database migrations.

Each time CAS Database models are updated, it is required to:

Generate a new migration:

.. code-block:: bash

    alembic -c casp/services/db/alembic.ini revision --autogenerate -m "{migration-message-goes-here}"

Apply migrations to the database:

.. code-block:: bash

    alembic -c casp/services/db/alembic.ini upgrade head
