Database Service
================

This module represents the codebase for connecting to the SQL database and performing operations on it. It is not a deployable service,
but rather a library that is used by other services.

.. note:: This service requires a running database cluster.

Setting up a Database Cluster
-----------------------------

Create a db cluster (`More info  <https://cloud.google.com/sdk/gcloud/reference/sql/instances/create>`_)

.. code-block:: bash

    DB_CLUSTER_NAME=cluster-example-name # Cluster name
    NUM_CPUS=4 # Number of CPUs
    MEMORY=15GB # Memory size
    STORAGE_SIZE=100GB # Disk storage size

    gcloud sql instances create $DB_CLUSTER_NAME \
    --database-version=POSTGRES_14 \
    --cpu=$NUM_CPUS \
    --memory=$MEMORY \
    --storage-size=$STORAGE_SIZE \
    --require-ssl \
    --region=us-central

Create a db user and password for the cluster (`More info  <https://cloud.google.com/sdk/gcloud/reference/sql/users/create>`_)

.. code-block:: bash

    DB_CLUSTER_NAME=cluster-example-name # Cluster name
    USERNAME=example-user # Username
    SECRET_PASSWORD='secretPassW0rd' # Password

    gcloud sql users create $USERNAME \
    --instance=$DB_CLUSTER_NAME \
    --password=$SECRET_PASSWORD

Create a database for the cluster (`More info  <https://cloud.google.com/sdk/gcloud/reference/sql/databases/create>`_)

.. code-block:: bash

    DB_NAME=example-db-name # Database name
    DB_CLUSTER_NAME=cluster-example-name # Cluster name

    gcloud sql databases create $DB_NAME \
    --instance=$DB_CLUSTER_NAME

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

Note: migration operations must be run from the cellarium-cloud/src directory.

Each time CAS Database models are updated, it is required to:

Generate a new migration:

.. code-block:: bash

    alembic -c casp/services/db/alembic.ini revision --autogenerate -m "{migration-message-goes-here}"

Note: You may need to modify the migration file in certain cases (e.g. backfilling data into new required columns)

Apply migrations to the database:

.. code-block:: bash

    alembic -c casp/services/db/alembic.ini upgrade head
