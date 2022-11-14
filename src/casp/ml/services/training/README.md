# Training Service Description
It is required to have a `src/casp/ml/services/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string
## Building Docker Image
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0-gpu
docker build -t $IMAGE_NAME -f ./src/casp/ml/services/deploy/Dockerfile.pytorch_cuda .
docker push $IMAGE_NAME
```
## Running Training Task via Cromshell
```
cd src/casp/ml/services/training/pca
cromshell submit train_command.wdl inputs_train_pca.json options.json
```
