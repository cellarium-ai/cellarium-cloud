# Inference Service Description
It is required to have a `src/casp/services/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string

## Building Docker Image
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0
docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch .
docker push $IMAGE_NAME
```
## Deploying Docker Image via Cloud Run
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0
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
--command python --args "casp/services/model_inference/server.py"
```