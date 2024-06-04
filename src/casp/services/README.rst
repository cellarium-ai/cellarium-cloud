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


It is required to have a ``src/settings/.env`` settings file.
The service must also be authenticated with Google credentials, either by setting:
- ``GOOGLE_APPLICATION_CREDENTIALS`` - path to the service account json credentials file
- Running in a Google environment (e.g. Cloud Run) with the service account attached to the service, in which case the service account is
used as the identity running the servce.

There are 2 CAS service types:

* Deployable solution: could be used in Cloud Run or deployed anywhere else
* Cromwell executable task

Both types are containerized solutions. There are 2 types of docker containers for CAS services:

#. Pytorch powered container
#. Pytorch + cuda powered container (needs a GPU to be executed)


Building Docker Images Locally
..............................


.. code-block:: bash

    IMAGE_NAME=docker-image.dev/example

    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch_cuda .
    docker push $IMAGE_NAME

Building Docker Images Remotely
...............................

This can be done by running the CAS Repository Docker Build workflow. This will build the docker image and push it to the CAS Repository.

This can be done via the github UI from:

https://github.com/cellarium-ai/cellarium-cloud/actions/workflows/docker-workflow.yml

Or by running with the Github client using the following command:

.. code-block:: bash

    GIT_REF=<branch, git commit, tag, etc.> # Git reference to build the image from
    IMAGE_TYPES=<standard|cuda|both> # Flavor of docker image to build: standard, cuda, both
    IMAGE_TAG=<docker tag> # Tag to apply to the docker image. If left empy, a short hash of the last git commit in $GIT_REF will be be used
    ADD_LATEST_TAG=<true|false> # If true, the image will be tagged as latest

    gh workflow run docker-workflow.yml --repo cellarium-ai/cellarium-cloud \
            --ref $GIT_REF \
            -f image-types=$IMAGE_TYPES \
            -f image-tag=$IMAGE_TAG \
            -f add-latest-tag=$ADD_LATEST_TAG

Deploying Docker Images Remotely
...............................

This can be done by running the CAS Repository Deploy workflow. This will deploy the docker image to Cloud Run.

This can be done via the github UI from:

https://github.com/cellarium-ai/cellarium-cloud/actions/workflows/deploy-workflow.yml

Or by running with the Github client using the following command:

.. code-block:: bash

    GIT_REF=<branch, git commit, tag, etc.> # Git reference to deploy the image from (helpful if modifying the deploy scripts or configuration files)
    IMAGE_TAG=<docker tag> # Tag to deploy the docker image from.
    SERVICE_ACCOUNT=sa-user@<project>.iam.gserviceaccount.com # Service account that will be running the service
    SECRET_REF=secret-name:latest # Reference to secret in the project's google secret manager as <secret name>:<version or latest> (note that the service account must have access to the secret)
    DB_CONNECTION=example-project:us-region-example:db-cluster-name # Cloud SQL connection name
    VPC_CONNECTOR=projects/<project>/locations/<region>/connectors/<connector> # VPC connector to use for the service
    DEPLOYMENT_PREFIX=example # Prefix to use for the deployment names
    DEPLOYMENT_FLAVOR=default # Flavor of the deployment to use. Must match a config file in the casp/services/deploy/configs directory

    gh workflow run deploy-workflow.yml repo cellarium-ai/cellarium-cloud \
            --ref $GIT_REF \
            -f image-tag=$IMAGE_TAG \
            -f service-account-email=$SERVICE_ACCOUNT \
            -f sql-instance=$DB_CONNECTION \
            -f vpc-connector=$VPC_CONNECTOR \
            -f config-secret=$SECRET_REF \
            -f deployment-prefix=$DEPLOYMENT_PREFIX \
            -f flavor=$DEPLOYMENT_FLAVOR


Service Architecture
....................


If it is a deployable service, it has to follow the following architecture:

#. ``main.py`` - entrypoint for the service (this is where FastAPI (or any other) application is initialized with configuration)
#. ``data_manager`` - module responsible for data access and communication with databases wherever the data's coming from
#. ``services`` - module responsible for domain logic (e.g. ModelInferenceService)
#. ``clients`` - module responsible for communication with other services (e.g. API communicates with model service)
#. ``schemas`` - module responsible for data validation using pydantic schemas
#. ``dependencies`` - module responsible for dependency injection (e.g. authentication)
