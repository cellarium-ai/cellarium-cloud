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

 - Python 3.10+

Developer Setup
~~~~~~~~~~~~~~~

To create a virtual python environment:

.. code-block:: shell

    python3 -mvenv python
    source python/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt
    pip install -r dev-requirements.txt
    pip install -e .


To run unit tests:

.. code-block:: shell

    tox -e unit

To lint:

.. code-block:: shell

    tox -e lint

To automatically fix formatting issues:

.. code-block:: shell

    tox -e format

Repository Versioning
---------------------
The repository uses `Semantic Versioning <https://semver.org/>`_ for versioning. For the versions available, see the
`tags on this repository <https://github.com/cellarium-ai/cellarium-cloud/tags>`_. Versioning is connected with the
repository branches. The `main` branch is the main production branch and has stable versions. The `development` branch
contains the latest development versions (such as release candidates (rc)). Other branches can be used for creating
alpha or beta versions.
