#!/bin/bash

# Initialize variables to default values
PROJECT_ID=""
REGION=""
INDEX_ENDPOINT_NAME=""
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
        --index_endpoint_name=*)
        INDEX_ENDPOINT_NAME="${arg#*=}"
        shift # Remove --index_endpoint_name= from processing
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

# This is a identifier and a display name YOU give for this deployed index (can be the same)
DEPLOYED_INDEX_ID="deployed_${INDEX_NAME}"
DISPLAY_NAME=$DEPLOYED_INDEX_ID

echo "Receiving endpoint_id..."
# Then we need the endpoint id with a little JQ magic
ENDPOINT_ID=$(gcloud ai index-endpoints list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_ENDPOINT_NAME\") | .name ")

echo "Receiving index_id..."
# and the id of the index to be deployed
INDEX_ID=$(gcloud ai indexes list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_NAME\") | .name ")

echo "Undeploying index from endpoint..."
# Undeploy Index from Endpoint
gcloud ai index-endpoints undeploy-index ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION} --deployed-index-id=${DEPLOYED_INDEX_ID}

echo "Removing index endpoint..."
gcloud ai index-endpoints delete ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION}

echo "Removing index..."
# Delete Index
gcloud ai indexes delete ${INDEX_ID} --project ${PROJECT_ID} --region ${REGION}
echo "Cleaned up!"