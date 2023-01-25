# Embedding Service Description
It is required to have a `src/casp/services/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string
## Building Docker Image
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0-gpu
docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch_cuda .
docker push $IMAGE_NAME
```
## Running Embedding Task via Cromshell
```
cd src/casp/services/embedding/pca
cromshell submit embedding_command.wdl inputs_embed_pca.json options.json
```
