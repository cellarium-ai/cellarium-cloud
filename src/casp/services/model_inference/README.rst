Model Inference Service
=======================
This service is responsible for running Cellarium Cloud model inference on a given dataset. It is deployed as a Cloud Run service.

Requirements
------------
- Python 3.10
- Pytorch 2.0.1
- Database connection
- `settings/.env` file with the secret variables


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
    MAX_INSTANCES=200 # Maximum number of instances to scale to
    MIN_INSTANCES=0 # Minimum number of instances to scale to. If 0, the service will have a "cold start"
    CONCURRENCY=10 # Maximum number of requests that can be served at the same time per instance
    SERVICE_ACCOUNT=sa-user@<project>.iam.gserviceaccount.com # Service account that will be running the service
    SECRET_REF=secret-name:latest # Reference to secret in the project's google secret manager as <secret name>:<version or latest> (note that the service account must have access to the secret)

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
    --service-account=$SERVICE_ACCOUNT \
    --set-secrets=/app/settings/.env=${SECRET_REF} \
    --command python --args "casp/services/model_inference/main.py" \
    --allow-unauthenticated \
    --ingress internal \
    --platform=managed

You can also deploy the services using the deploy-workflow.yml GitHub action.

Accessing the Service
---------------------
The service is behind a firewall and not exposed to the public internet.  It can only be accessed by the internal or if you create 
a tunnel through the bastion host. To create a tunnel:

log into the gcloud cli with the command:
.. code-block:: bash

    gcloud auth login

then run the following command:

.. code-block:: bash

    gcloud compute ssh --zone "us-central1-a" "bastion" --project "dsp-cell-annotation-service" --ssh-flag="-D 9090" --ssh-flag="-N"

To use the proxy, you can configure your browser to use the SOCKS proxy at localhost:9090 (or whatever port you specified in the command above).

When you are done, you can close the tunnel by stopping the ssh command in the terminal.
