CAS Backend
===========

CAS Backend is the backend that powers CAS, the Cell Annotation Service, including the APIs used by clients and the
internal administration tooling that supports the service.
The repository is organized as a Python package rooted at ``cellarium.cas_backend`` and is documented using a
README-first structure: major modules keep their own ``README.rst`` files next to the code, and Sphinx assembles
those module READMEs into the published documentation at
`cellarium-cloud.readthedocs.io <https://cellarium-cloud.readthedocs.io>`_.

Repository Structure
--------------------

Top-level layout:

.. code-block:: text

    cellarium/cas_backend/
      apps/          deployable services and app-specific code
      core/          shared runtime, configuration, auth, data access, and db code
    deploy/          Dockerfiles, Cloud Run flavor configs, and deployment scripts
    docs/source/     Sphinx docs, including wrappers around module READMEs
    settings/        local-only environment files such as sample_env and .env
    tests/           unit and other automated tests

Documentation Contract
----------------------

This repository uses two documentation layers:

- Module ``README.rst`` files are the canonical source for module purpose, structure, entrypoints, dependencies, and
  local development notes.
- Sphinx pages in ``docs/source/`` own cross-cutting topics such as local development, CI/CD, secrets, deployment,
  and navigation across modules.

If you change a major module boundary or developer workflow, update the nearby module README first and then adjust any
cross-cutting Sphinx pages that reference it.

Prerequisites
-------------

- Python 3.12+
- Poetry

Developer Setup
---------------

Install Poetry:

.. code-block:: shell

    make setup

Install project dependencies, including development and test groups:

.. code-block:: shell

    make install

Create a local environment file from the sample and fill in the required secrets:

.. code-block:: shell

    cp settings/sample_env settings/.env

Common Commands
---------------

.. code-block:: shell

    make help
    make test
    make lint
    make format
    make requirements
    make docs

You can also run the underlying Poetry commands directly:

.. code-block:: shell

    poetry install --with dev,test
    poetry run pytest tests/unit
    poetry run ruff check cellarium tests
    poetry run ruff format cellarium tests

Documentation
-------------

Build the docs locally with:

.. code-block:: shell

    poetry install --with docs
    poetry run sphinx-build -W -b html docs/source docs/build/html

Start with these module entry points when navigating the codebase:

- ``cellarium/cas_backend/apps/README.rst``
- ``cellarium/cas_backend/core/README.rst``
- ``deploy/cloudrun/README.rst``

Versioning
----------

The repository uses `Semantic Versioning <https://semver.org/>`_. Stable releases are published from ``main``.
