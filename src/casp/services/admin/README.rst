Admin Service
=============

Admin service is a service that is used to manage the database and internal state of the application. Currently, all deployed models and indexes are stored in the database.
The service is built on top of `Flask <https://flask.palletsprojects.com/>`_ and uses a database connection to perform operations on the database.

Requirements
------------
- Python 3.10
- Database connection
- `settings/.env` file with secret variables


Building Docker Image
---------------------

.. code-block:: bash

    IMAGE_NAME=docker-image.dev/example # Name of the docker image
    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
    docker push $IMAGE_NAME

Deploying Docker Image via Cloud Run
------------------------------------

To deploy the Docker image using Cloud Run run (see `Cloud Run Documentation <https://cloud.google.com/sdk/gcloud/reference/run/deploy>`_ for more information)

.. code-block:: bash

    SERVICE_NAME=cas-admin-service # Name of the service
    PROJECT_ID=example-project # Name of the project
    IMAGE_NAME=docker-image.dev/example # Name of the docker image # Name of the docker image
    REGION=us-central1 # Region where the service will be deployed
    PORT=8000 # Port on which the service will be running (matches the port in the flask app)
    DB_CONNECTION=example-project:us-region-example:db-cluster-name # Cloud SQL connection name
    SERVICE_ACCOUNT=sa-user@<project>.iam.gserviceaccount.com # Service account that will be running the service
    SECRET_REF=secret-name:latest # Reference to secret in the project's google secret manager as <secret name>:<version or latest> (note that the service account must have access to the secret)

    gcloud run deploy $SERVICE_NAME \
    --project=$PROJECT_ID \
    --image=$IMAGE_NAME \
    --region=$REGION \
    --port=$PORT \
    --add-cloudsql-instances=$DB_CONNECTION \
    --service-account=$SERVICE_ACCOUNT \
    --set-secrets=/app/settings/.env=${SECRET_REF} \
    --command=casp/services/admin/entrypoint.sh \
    --platform managed \
    --allow-unauthenticated

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
