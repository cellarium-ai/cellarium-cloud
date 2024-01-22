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

.. code-block:: bash

    IMAGE_NAME=us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.4.1
    docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
    docker push $IMAGE_NAME

Deploying Docker Image via Cloud Run
------------------------------------

.. code-block:: bash

    IMAGE_NAME=us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.4.1
    PROJECT_ID=dsp-cell-annotation-service

    gcloud run deploy cas-model \
    --project=$PROJECT_ID \
    --image=$IMAGE_NAME \
    --memory=4Gi \
    --cpu=1 \
    --region=us-central1 \
    --platform=managed \
    --port=8000 \
    --allow-unauthenticated \
    --add-cloudsql-instances=dsp-cell-annotation-service:us-central1:cas-db-cluster-2 \
    --timeout=1100 \
    --max-instances=500 \
    --min-instances=0 \
    --concurrency=10 \
    --command python --args "casp/services/model_inference/main.py"
