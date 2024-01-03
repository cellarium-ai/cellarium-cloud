# API Internal Service Description
This service is responsible for handling all internal requests from other services. It is not exposed to the outside 
world. It can be used to handle requests from the pipelines, or from the other services.

It is required to have a `src/casp/services/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string

## Building Docker Image
```
IMAGE_NAME=us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.4.2.dev
docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
docker push $IMAGE_NAME
```
## Deploying Docker Image via Cloud Run
```
IMAGE_NAME=us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.4.2.dev
PROJECT_ID=dsp-cell-annotation-service

gcloud run deploy cas-api-internal \
--project=$PROJECT_ID \
--image=$IMAGE_NAME \
--memory=1Gi \
--cpu=2 \
--region=us-central1 \
--platform=managed \
--port=8000 \
--add-cloudsql-instances=dsp-cell-annotation-service:us-central1:cas-db-cluster-dev \
--allow-unauthenticated \
--vpc-connector=cas-ai-matching \
--ingress=internal \
--max-instances=10 \
--min-instances=0 \
--concurrency=40 \
--command python --args "casp/services/api_internal/main.py"
```