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

PROJECT_NUMBER=`gcloud projects describe $PROJECT_ID | grep projectNumber | cut -d"'" -f2`

# This is a indentifier and a display name YOU give for this deployed index (can be the same)
DEPLOYED_INDEX_ID="deployed_${INDEX_NAME}"
DISPLAY_NAME=$DEPLOYED_INDEX_ID

# Then we need the endpoint id with a little JQ magic
ENDPOINT_ID=$(gcloud ai index-endpoints list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_ENDPOINT_NAME\") | .name ")

# and the id of the index to be deployed
INDEX_ID=$(gcloud ai indexes list --region $REGION --project $PROJECT_ID --format json | jq -r ".[] | select (.displayName == \"$INDEX_NAME\") | .name ")

gcloud ai index-endpoints deploy-index $ENDPOINT_ID \
  --deployed-index-id=$DEPLOYED_INDEX_ID \
  --display-name=$DISPLAY_NAME \
  --index=$INDEX_ID \
  --min-replica-count 1 \
  --max-replica-count 2 \
  --region=$REGION \
  --project=$PROJECT_ID
