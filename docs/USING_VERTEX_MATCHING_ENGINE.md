# Vertex AI Matchine Enging

## Summary - The good, the bad and the ugly

### Good
- Autoscaling (min/max nodes) is a nice feature (assuming it works)
- Has ability to restrict matches based on user-defined labels on the vectors
- Potentially supports incremental updates

### Bad
- index creation, even for 500 cells, took quite a long time (like 30 minutes!)

### Ugly
- UI is missing a lot of features
- network configuration is needlessly complicated
- some calls use project_id and some use project_number, some use fully-qualified paths other short names

## Configuration, Creating and Deploying an Index

Using Vertex AI Matching Engine requires several steps before matching

1. Network Configuration (one time)
2. Create Endpoint (one time)
3. Create Index
4. Deploy Index

### Constants

```
# Googe project, project number and region to host Vertex AI
export PROJECT_ID="dsp-cell-annotation-service"
export PROJECT_NUMBER=`gcloud projects describe $PROJECT_ID | grep projectNumber | cut -d"'" -f2`
export REGION="us-central1"

# Bucket containing CSV/AVRO of vectors to be searched
export BUCKET_URI="gs://dsp-cell-annotation-service/kcibul/matching/"
export DIMENSIONS=75
export APPROX_NEIGHBORS_COUNT=100

# Constants, not necessary to change
export VPC_NETWORK="ai-matching"
export PEERING_RANGE_NAME="ann-haystack-range"
export INDEX_ENDPOINT_NAME="casp_index_endpoint"
export INDEX_NAME="casp_index_v1"
```

### Network Configuration

In order to use Vertex AI Matching Engine, a VPC Network with Peering enabled must be created.  

#### Create the VPC Network
```
gcloud compute networks create ${VPC_NETWORK} --bgp-routing-mode=regional --subnet-mode=auto --project=${PROJECT_ID}
```

#### Add necessary firewall rules
```
gcloud compute firewall-rules create ${VPC_NETWORK}-allow-icmp --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow icmp

gcloud compute firewall-rules create ${VPC_NETWORK}-allow-internal --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow all --source-ranges 10.128.0.0/9

gcloud compute firewall-rules create ${VPC_NETWORK}-allow-rdp --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow tcp:3389

gcloud compute firewall-rules create ${VPC_NETWORK}-allow-ssh --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow tcp:22
```

#### Reserve IP range
```
gcloud compute addresses create ${PEERING_RANGE_NAME} --global --prefix-length=16 --network=${VPC_NETWORK} --purpose=VPC_PEERING --project=${PROJECT_ID}
```

#### Set up peering with service networking

Note: Your account must have the "Compute Network Admin" role to run the following.

```
gcloud services vpc-peerings connect --service=servicenetworking.googleapis.com --network=${VPC_NETWORK} --ranges=${PEERING_RANGE_NAME} --project=${PROJECT_ID}
```

### Create Index Endpoint (to serve the index)

This step takes several minutes to complete

```
gcloud ai index-endpoints create --display-name ${INDEX_ENDPOINT_NAME} --network projects/${PROJECT_NUMBER}/global/networks/${VPC_NETWORK} --region ${REGION} --project $PROJECT_ID
```

### Create Index

Create the actual index, this takes a long time! (~30 minutes even for a small dataset)

```
# save configuration to a local file
export LOCAL_PATH_TO_METADATA_FILE=/tmp/metadata.json
cat << EOF > ${LOCAL_PATH_TO_METADATA_FILE}
{
  "contentsDeltaUri": "${BUCKET_URI}",
  "config": {
    "dimensions": ${DIMENSIONS},
    "approximateNeighborsCount": ${APPROX_NEIGHBORS_COUNT},
    "distanceMeasureType": "DOT_PRODUCT_DISTANCE",
    "algorithm_config": {
      "treeAhConfig": {
      }
    }
  }
}
EOF

gcloud ai indexes create \
  --metadata-file=${LOCAL_PATH_TO_METADATA_FILE} \
  --display-name=${INDEX_NAME} \
  --project=${PROJECT_ID} \
  --region=${REGION}

```

This is an async operation, you'll have to poll for success (the command is given by the create command above)

For example, 

```
gcloud ai operations describe 2843220864793575424 --index=7139735929568100352 --region us-central1 --project=dsp-cell-annotation-service
```

### Deploy Index

Deploy the index to the endpoint so it can be searched.  Several non-intuitive IDs are required to run this step

```
# This is a indentifier and a display name YOU give for this deployed index (can be the same)
export DEPLOYED_INDEX_ID="deployed_${INDEX_NAME}"
export DISPLAY_NAME=$DEPLOYED_INDEX_ID

# Then we need the endpoint id with a little JQ magic
export ENDPOINT_ID=$(gcloud ai index-endpoints list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_ENDPOINT_NAME\") | .name ")

# and the id of the index to be deployed
export INDEX_ID=$(gcloud ai indexes list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_NAME\") | .name ")

gcloud ai index-endpoints deploy-index $ENDPOINT_ID \
  --deployed-index-id=$DEPLOYED_INDEX_ID \
  --display-name=$DISPLAY_NAME \
  --index=$INDEX_ID \
  --min-replica-count 2 \
  --max-replica-count 2 \

```

This is an async operation, you'll have to poll for success (the command is given by the create command above)

For example, 
```
gcloud ai operations describe 1574402038526115840 --index-endpoint=82032363525111808 --project $PROJECT_ID --region $REGION
```

## Search!

Searching can only be performed from compute on the same network that was configured about with the proper peering settings.  The easiest way to do this is to create a Notebook instance and under the Networking configuration choose the VPC network created in the above steps (i.e. `ai-matching` in this example).

The DIMENSIONS, ENDPOINT_ID and DEPLOYED_INDEX_ID variables should have the value from above

Then from that notebook VM:

```
from google.cloud import aiplatform
import numpy as np

DIMENSIONS=75
ENDPOINT_ID="projects/350868384795/locations/us-central1/indexEndpoints/82032363525111808"
DEPLOYED_INDEX_ID="deployed_casp_index_v1"

# locate the endpoint
ep = aiplatform.MatchingEngineIndexEndpoint(index_endpoint_name=ENDPOINT_ID)

# generate a random vector to search with
emb1 = np.random.randn(75)

# perform the query
response = index_endpoint.match(deployed_index_id=DEPLOYED_INDEX_ID, queries=[emb1], num_neighbors=25)

# response is an array of results where each result is an array of MatchNeighbor objects
for result in response:
    for match in result:
        print(f"ID:{match.id} DISTANCE:{match.distance}")
```

## Evaluating Performance

Aspects to consider:

1. Throughput (overall matches per second)
2. Latency (response time per request)
3. Scalability (with respect to index size)
4. Accuracy
5. Cost

TBD

## Cleaning Up (excluding the network setup)

If you want to remove everything, just go in the opposite order from the above

```
# Undeploy Index from Endpoint
gcloud ai index-endpoints undeploy-index ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION} --deployed-index-id=${DEPLOYED_INDEX_ID}

# Delete Endpoint
gcloud ai index-endpoints delete ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION}

# Delete Index
gcloud ai indexes delete ${INDEX_ID} --project ${PROJECT_ID} --region ${REGION}
```