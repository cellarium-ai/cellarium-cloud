Model Inference Service
=======================
This service is responsible for running Cellarium Cloud model inference on a given dataset. It is deployed as a Cloud Run service.

Requirements
------------
- Python 3.10
- Pytorch 2.0.1
- Database connection
- `.env` file with the secret variables


Building Docker Image
---------------------
To build the Docker image locally and push it to the registry, run the following commands:

.. code-block:: bash

    IMAGE_NAME=docker-image.dev/example # Name of the docker image
    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
    docker push $IMAGE_NAME

Deploying Docker Image via Cloud Run
------------------------------------

To deploy the Docker image using Cloud Run run (see `Cloud Run Documentation <https://cloud.google.com/sdk/gcloud/reference/run/deploy>`_ for more information)

.. code-block:: bash

    SERVICE_NAME=cellarium-cloud-model-inference # Name of the service
    IMAGE_NAME=docker-image.dev/example # Name of the docker image
    PROJECT_ID=example-project # GCP Project ID
    NUM_CPUS=4 # Number of CPUs per deployed instance
    MEMORY=16Gi # Memory per deployed instance; Can't be less than 256Mi and more than 4Gi per one CPU core
    REGION=us-central1 # Region where the service will be deployed
    PLATFORM=managed # Target platform to run the service. Choices: managed, gke, kubernetes
    PORT=8000 # Port which the running image will listen to (matches the FastAPI port).
    DB_CONNECTION=example-project:us-region-example:db-cluster-name # Cloud SQL connection name
    TIMEOUT=1100 # Request timeout in seconds
    MAX_INSTANCES=500 # Maximum number of instances to scale to
    MIN_INSTANCES=0 # Minimum number of instances to scale to. If 0, the service will have a "cold start"
    CONCURRENCY=10 # Maximum number of requests that can be served at the same time per instance

    gcloud run deploy $SERVICE_NAME \
    --project=$PROJECT_ID \
    --image=$IMAGE_NAME \
    --memory=$MEMORY \
    --cpu=$NUM_CPUS \
    --region=$REGION \
    --port=$PORT \
    --add-cloudsql-instances=$DB_CONNECTION \
    --timeout=$TIMEOUT \
    --max-instances=$MAX_INSTANCES \
    --min-instances=$MIN_INSTANCES \
    --concurrency=$CONCURRENCY \
    --command python --args "casp/services/model_inference/main.py" \
    --allow-unauthenticated \
    --platform=managed
