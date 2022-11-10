# Training Service Description
It is required to have a `src/casp/ml/services/training/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string
## Building Docker Image
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_training-gpu:1.0
docker build -t $IMAGE_NAME -f ./src/casp/ml/services/training/Dockerfile.training_service .
docker push $IMAGE_NAME
```
## Running Training Task via Cromshell
```
cromshell submit \
./src/casp/ml/services/training/train_command.wdl \
./src/casp/ml/services/training/inputs.json \ 
./src/casp/ml/services/training/options.json
```
