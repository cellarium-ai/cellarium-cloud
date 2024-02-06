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
- `.env` file with the secret variables

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

To build and push the Docker image:

.. code-block:: bash

   IMAGE_NAME=us-docker.pkg.dev/dsp-cell-annotation-service/cas/cas-api:0.1
   docker build -t $IMAGE_NAME -f Dockerfile.api .
   docker push $IMAGE_NAME

Set up VPC
----------
.. note::

    This has to be executed only once.

In order for Cloud Run services to access the ai-matching network, which is necessary to access the Vertex AI Matching Engine, Serverless VPC access must be enabled. This comes with a small cost to run the e2-micro VMs that do the network bridging.

.. code-block:: bash

   gcloud compute networks vpc-access connectors create cas-ai-matching \
   --project dsp-cell-annotation-service \
   --region=us-central1 \
   --network=ai-matching \
   --range=10.8.0.0/28 \
   --min-instances=2 \
   --max-instances=10 \
   --machine-type=e2-micro

Deploying Docker Image via Cloud Run
------------------------------------

To deploy the Docker image using Cloud Run:

.. code-block:: bash

   IMAGE_NAME=us-docker.pkg.dev/dsp-cell-annotation-service/cas/cas-api:0.1
   PROJECT_ID=dsp-cell-annotation-service

   gcloud run deploy cas-api \
   --project $PROJECT_ID \
   --image $IMAGE_NAME \
   --cpu=1 \
   --memory=4Gi \
   --region=us-central1 \
   --platform=managed \
   --port=8000 \
   --allow-unauthenticated \
   --vpc-connector=cas-ai-matching \
   --add-cloudsql-instances=dsp-cell-annotation-service:us-central1:cas-db-cluster-2 \
   --timeout=1100 \
   --max-instances=500 \
   --min-instances=0 \
   --concurrency=20 \
   --command=python --args="casp/services/api/main.py"

Test your deployment with:

.. code-block:: bash

   BASE_URL="https://cas-api-vi7nxpvk7a-uc.a.run.app:8000"
   curl -X POST -H "Accept: application/json" -F "json=\"gimme-som-data\";type=application/json" -F "myfile=@local_1000.h5ad" "$BASE_URL/annotate" -o results.json
