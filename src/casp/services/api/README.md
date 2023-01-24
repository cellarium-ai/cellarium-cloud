## Running Locally

Note: VM Must be on the network as Vertex Matching Engine (see docs/USING_VERTEX_MATCHING_ENGINE.md)

```
pip install -r requirements.txt
python src/casp/services/api/server.py
```

## Building Docker Image

```
IMAGE_NAME=us-docker.pkg.dev/dsp-cell-annotation-service/cas/cas-api:0.1
docker build -t $IMAGE_NAME -f Dockerfile.api .
docker push $IMAGE_NAME
```
## Set up VPC 
In order for Cloud Run services to access the ai-matching network, which is necessary to access the Vertex AI Matching Engine, Serverless VPC access must be enabled.  This comes with a small cost to run the e2-micro VMs that do the network bridging.

```
gcloud compute networks vpc-access connectors create cas-ai-matching \
--project dsp-cell-annotation-service \
--region=us-central1 \
--network=ai-matching \
--range=10.8.0.0/28 \
--min-instances=2 \
--max-instances=10 \
--machine-type=e2-micro
```

## Deploying Docker Image via Cloud Run
```
IMAGE_NAME=us-docker.pkg.dev/dsp-cell-annotation-service/cas/cas-api:0.1
PROJECT_ID=dsp-cell-annotation-service

gcloud run deploy cas-api-test-2 \
--project $PROJECT_ID \
--image $IMAGE_NAME \
--memory 4Gi \
--region us-central1 \
--platform managed \
--port 8000 \
--allow-unauthenticated \
--vpc-connector cas-ai-matching \
--command python --args "casp/services/api/server.py"


BASE_URL="https://cas-api-vi7nxpvk7a-uc.a.run.app:8000"
curl -X POST -H "Accept: application/json" -F "json=\"gimme-som-data\";type=application/json" -F "myfile=@local_1000.h5ad" "$BASE_URL/annotate" -o results.json
