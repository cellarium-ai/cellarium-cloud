Running Locally
===============

This page covers the cross-cutting local setup for CAS Backend. Module-specific details stay in the module READMEs.

Prerequisites
-------------

- Python 3.12
- Poetry
- A local or reachable PostgreSQL instance
- A populated ``settings/.env`` file at the repository root

Initial Setup
-------------

.. code-block:: shell

    make setup
    make install
    cp settings/sample_env settings/.env

Fill in the values required for your environment before starting services.

Running the Compute Service
---------------------------

.. code-block:: shell

    poetry run python -m cellarium.cas_backend.apps.compute.main

OpenAPI docs are available at `http://localhost:8000/api/docs <http://localhost:8000/api/docs>`_.

Running the Admin Service
-------------------------

.. code-block:: shell

    poetry run python -m cellarium.cas_backend.apps.admin.server

Testing and Linting
-------------------

.. code-block:: shell

    make test
    make lint
