API Service
===========

The API service is a REST API that provides access to the Cellarium Backend. It is built using FastAPI and deployed
using Cloud Run. It is designed to be deployed in the same VPC as the Vertex Matching Engine as this is the only way how
the API can access the Vertex Matching Engine. It is also required to have a cloud SQL connection to the database.

API is the only service that communicates with the outside world. Its methods require authentication and authorization.

Requirements
------------
- Python 3.10
- Database connection
- Vertex Matching Engine connection
- `settings/.env` file with the secret variables

Running Locally
---------------

.. note::

   VM Must be on the network as Vertex Matching Engine (see :doc:`Vertex Matching Engine Documentation <../vertex_ai_matching_engine>`)

To install dependencies and run the server locally:

.. code-block:: bash

   pip install -r requirements.txt
   python src/casp/services/api/main.py

Building Docker Image
---------------------

To build the Docker image locally and push it to the registry, run the following commands:

.. code-block:: bash

   IMAGE_NAME=docker-image.dev/example
   docker build -t $IMAGE_NAME -f Dockerfile.api .
   docker push $IMAGE_NAME

Set up VPC
----------
.. note::

    This has to be executed only once.

In order for Cloud Run services to access the VPC network, which is necessary to access the Vertex AI Matching Engine, Serverless VPC access must be enabled. This comes with a small cost to run the e2-micro VMs that do the network bridging.

To create a VPC connector run the following command (`More info <https://cloud.google.com/vpc/docs/configure-serverless-vpc-access>`_):

.. code-block:: bash

    VPC_CONNECTOR_NAME=example-vpc-connector-name # Name of the VPC connector
    PROJECT_ID=example-project # GCP Project ID
    REGION=us-central1 # Region where the VPC connector will be deployed
    NETWORK_NAME=ai-matching # Name of the network to connect to
    IP_RANGE=10.8.0.0/28 # IP range for the VPC connector
    MIN_INSTANCES=2 # Minimum number of instances to run
    MAX_INSTANCES=10 # Maximum number of instances to run

    gcloud compute networks vpc-access connectors create $VPC_CONNECTOR_NAME \
    --project=$PROJECT_ID \
    --region=$REGION \
    --network=$NETWORK_NAME \
    --range=$IP_RANGE \
    --min-instances=$MIN_INSTANCES \
    --max-instances=$MAX_INSTANCES \
    --machine-type=e2-micro

Deploying Docker Image via Cloud Run
------------------------------------

To deploy the Docker image using Cloud Run run (see `Cloud Run Documentation <https://cloud.google.com/sdk/gcloud/reference/run/deploy>`_ for more information)

.. code-block:: bash

    SERVICE_NAME=cellarium-cloud-api # Name of the service
    IMAGE_NAME=docker-image.dev/example # Name of the docker image
    PROJECT_ID=example-project # GCP Project ID
    NUM_CPUS=1 # Number of CPUs per deployed instance
    MEMORY=4Gi # Memory per deployed instance; Can't be less than 256Mi and more than 4Gi per one CPU core
    REGION=us-central1 # Region where the service will be deployed
    PLATFORM=managed # Target platform to run the service. Choices: managed, gke, kubernetes
    PORT=8000 # Port which the running image will listen to (matches the FastAPI port).
    VPC_CONNECTOR=example-vpc-connector-name # Name of the VPC connector (matches the name of the VPC connector for the Vertex Matching Engine)
    DB_CONNECTION=example-project:us-region-example:db-cluster-name # Cloud SQL connection name
    TIMEOUT=1100 # Request timeout in seconds
    MAX_INSTANCES=200 # Maximum number of instances to scale to
    MIN_INSTANCES=0 # Minimum number of instances to scale to. If 0, the service will have a "cold start"
    CONCURRENCY=20 # Maximum number of requests that can be served at the same time per instance
    SERVICE_ACCOUNT=sa-user@<project>.iam.gserviceaccount.com # Service account that will be running the service
    SECRET_REF=secret-name:latest # Reference to secret in the project's google secret manager as <secret name>:<version or latest> (note that the service account must have access to the secret)

    gcloud run deploy $SERVICE_NAME \
    --project=$PROJECT_ID \
    --image=$IMAGE_NAME \
    --cpu=$NUM_CPUS \
    --memory=$MEMORY \
    --region=$REGION \
    --port=$PORT \
    --timeout=$TIMEOUT \
    --max-instances=$MAX_INSTANCES \
    --min-instances=$MIN_INSTANCES \
    --concurrency=$CONCURRENCY \
    --add-cloudsql-instances=$DB_CONNECTION \
    --vpc-connector=$VPC_CONNECTOR \
    --platform=$PLATFORM \
    --service-account=$SERVICE_ACCOUNT \
    --set-secrets=/app/settings/.env=${SECRET_REF} \
    --command=python --args="casp/services/api/main.py" \
    --allow-unauthenticated

Test your deployment with:

.. code-block:: bash

   BASE_URL="https://the-url-of-the-deployed-service.dev"
   curl -X POST -H "Accept: application/json" -F "json=\"gimme-som-data\";type=application/json" -F "myfile=@local_1000.h5ad" "$BASE_URL/annotate" -o results.json
