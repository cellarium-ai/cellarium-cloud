# CAS Services
It is required to have a `src/casp/services/.env` with `GOOGLE_SERVICE_ACCOUNT_CREDENTIALS` variable which includes service account json credentials dumped as oneline string to use any services.
There are 2 CAS service types:
* Deployable solution: could be used in Cloud Run or deployed anywhere else
* Cromwell executable task

Both types are containerized solutions. \
There are 2 types of docker containers for CAS services:
1. Pytorch powered container
2. Pytorch + cuda powered container (needs a GPU to be executed)

## Building Docker Images
### 
```
IMAGE_NAME=us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_service:1.0-gpu
docker build -t $IMAGE_NAME -f ./src/casp/services/deploy/Dockerfile.pytorch_cuda .
docker push $IMAGE_NAME
```

## Service Architecture
If it is a deployable service, it has to follow the following architecture:
1. `main.py` - entrypoint for the service (this is where FastAPI (or any other) application is initialized with configuration)
2. `data_manager` - module responsible for data access and communication with databases wherever the data's coming from
3. `services` - module responsible for domain logic (e.g. ModelInferenceService)
4. `clients` - module responsible for communication with other services (e.g. API communicates with model service)
5. `schemas` - module responsible for data validation using pydantic schemas
6. `dependencies` - module responsible for dependency injection (e.g. authentication)