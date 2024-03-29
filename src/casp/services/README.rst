Cellarium Cloud Services
========================

This section provides an overview of the various services within the Cellarium Cloud platform.


Below is the list of individual services:

.. toctree::
   :maxdepth: 1

   api_service
   model_inference_service
   admin_service
   db_service

General Description
-------------------


It is required to have a ``src/casp/services/.env`` with ``GOOGLE_SERVICE_ACCOUNT_CREDENTIALS`` variable which includes service account json credentials dumped as a one-line string to use any services.
There are 2 CAS service types:

* Deployable solution: could be used in Cloud Run or deployed anywhere else
* Cromwell executable task

Both types are containerized solutions. There are 2 types of docker containers for CAS services:

#. Pytorch powered container
#. Pytorch + cuda powered container (needs a GPU to be executed)


Building Docker Images
......................


.. code-block:: bash

    IMAGE_NAME=docker-image.dev/example

    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch_cuda .
    docker push $IMAGE_NAME

Service Architecture
....................


If it is a deployable service, it has to follow the following architecture:

#. ``main.py`` - entrypoint for the service (this is where FastAPI (or any other) application is initialized with configuration)
#. ``data_manager`` - module responsible for data access and communication with databases wherever the data's coming from
#. ``services`` - module responsible for domain logic (e.g. ModelInferenceService)
#. ``clients`` - module responsible for communication with other services (e.g. API communicates with model service)
#. ``schemas`` - module responsible for data validation using pydantic schemas
#. ``dependencies`` - module responsible for dependency injection (e.g. authentication)
