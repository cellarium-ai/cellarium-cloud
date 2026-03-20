#!/bin/bash

# Initialize variables to default values
PROJECT_ID=""
REGION=""
BUCKET_URI=""
DIMENSIONS=""
APPROX_NEIGHBORS_COUNT=""
INDEX_NAME=""

# Loop through all the passed arguments
for arg in "$@"
do
    case $arg in
        --project_id=*)
        PROJECT_ID="${arg#*=}"
        shift # Remove --project_id= from processing
        ;;
        --region=*)
        REGION="${arg#*=}"
        shift # Remove --region= from processing
        ;;
        --bucket_uri=*)
        BUCKET_URI="${arg#*=}"
        shift # Remove --bucket_uri= from processing
        ;;
        --dimensions=*)
        DIMENSIONS="${arg#*=}"
        shift # Remove --dimensions= from processing
        ;;
        --approx_neighbors_count=*)
        APPROX_NEIGHBORS_COUNT="${arg#*=}"
        shift # Remove --approx_neighbors_count= from processing
        ;;
        --index_name=*)
        INDEX_NAME="${arg#*=}"
        shift # Remove --index_name= from processing
        ;;
        *)
        # Unknown option
        echo "Unknown argument $arg"
        exit 1
        ;;
    esac
done

echo "Receiving project number..."
PROJECT_NUMBER=`gcloud projects describe $PROJECT_ID | grep projectNumber | cut -d"'" -f2`

echo "Creating Index from embeddings provided via 'bucket_uri'..."
# save configuration to a local file
LOCAL_PATH_TO_METADATA_FILE=/tmp/metadata.json
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

echo "Job Submitted. You can log into Google Vertex AI dashboard for details. It should take ~1 hr to create an index."