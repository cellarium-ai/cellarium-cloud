Cellarium Cloud
===============

**What's Cellarium Cloud?** Cellarium Cloud contains the service code that supports the Cellarium Cloud platform by providing APIs for the Cellarium Client tool.
The APIs can be access with direct REST access or, preferably, via the cellarium-cas Python client.

Repository Structure
--------------------
Each module of the repository contains its own `README.rst` file that explains the purpose of the module and how to run
it. These files are united together in a documentation website that is available at
`Cellarium Cloud Documentation <https://cellarium-cloud.readthedocs.io>`_.


Prerequisites / Installation
----------------------------

 - Python 3.12+
 - Poetry (Python dependency management tool)

Developer Setup
~~~~~~~~~~~~~~~

To set up the development environment, first install Poetry:

.. code-block:: shell

    make setup

Then install all dependencies (including dev and test dependencies):

.. code-block:: shell

    make install

Available Make Commands
~~~~~~~~~~~~~~~~~~~~~~~

The project uses Make commands to streamline common development tasks:

**Setup and Installation:**

.. code-block:: shell

    make setup    # Install Poetry
    make install  # Install all dependencies with Poetry

**Testing:**

.. code-block:: shell

    make test     # Run unit tests with pytest

**Code Quality:**

.. code-block:: shell

    make lint     # Run ruff linting (check only)
    make format   # Auto-format code with ruff

**Dependency Management:**

.. code-block:: shell

    make requirements        # Export production dependencies to deploy/requirements.txt.lock
    make docker-requirements # Export all dependencies for Docker builds

**Cleanup:**

.. code-block:: shell

    make clean    # Remove generated files (coverage, cache, etc.)

For a full list of available commands, run:

.. code-block:: shell

    make help

Using Poetry Directly
~~~~~~~~~~~~~~~~~~~~~

You can also use Poetry commands directly:

.. code-block:: shell

    poetry install --with dev,test           # Install dependencies
    poetry run pytest tests/unit             # Run tests
    poetry run ruff check cellarium tests    # Lint code
    poetry run ruff format cellarium tests   # Format code

Repository Versioning
---------------------
The repository uses `Semantic Versioning <https://semver.org/>`_ for versioning. For the versions available, see the
`tags on this repository <https://github.com/cellarium-ai/cellarium-cloud/tags>`_. Versioning is connected with the
repository branches. The `main` branch is the main production branch and has stable versions. The `development` branch
contains the latest development versions (such as release candidates (rc)). Other branches can be used for creating
alpha or beta versions.
