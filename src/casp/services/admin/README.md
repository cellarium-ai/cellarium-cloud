# Admin Service Description
## Depends on `db` module
It is required to have a `src/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string

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

gcloud run deploy casp-admin-service \
--project $PROJECT_ID \
--image $IMAGE_NAME \
--region us-central1 \
--platform managed \
--port 8000 \
--allow-unauthenticated \
--add-cloudsql-instances=dsp-cell-annotation-service:us-central1:cas-db-cluster \
--command casp/services/admin/entrypoint.sh
```