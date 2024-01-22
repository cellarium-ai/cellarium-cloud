Admin Service
=============

Admin service is a service that is used to manage the database and internal state of the application. Currently, all deployed models and indexes are stored in the database.
The service is built on top of `Flask <https://flask.palletsprojects.com/>`_ and uses a database connection to perform operations on the database.

Requirements
------------
- Python 3.10
- Database connection
- `.env` file with secret variables


Building Docker Image
---------------------

.. code-block:: bash

    IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0
    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
    docker push $IMAGE_NAME

Deploying Docker Image via Cloud Run
------------------------------------

.. code-block:: bash

    IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0
    ZONE=us-central1
    DB_CONNECTION=dsp-cell-annotation-service:us-central1:cas-db-cluster
    PROJECT_ID=dsp-cell-annotation-service

    gcloud run deploy casp-admin-service \
    --project $PROJECT_ID \
    --image $IMAGE_NAME \
    --region $ZONE \
    --platform managed \
    --port 8000 \
    --allow-unauthenticated \
    --add-cloudsql-instances=$DB_CONNECTION \
    --command casp/services/admin/entrypoint.sh
