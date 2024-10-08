Vertex AI Matching Engine
=========================

Vertex AI is a Google Cloud service that provides a variety of machine learning tools. One of these tools is the Matching Engine, 
which is a vector similarity search service. The Matching Engine is used in CAS to search for cells that are in close proximity in a vector space to the user's submitted cells.

Summary - The good, the bad and the ugly
----------------------------------------

Good
^^^^
- Autoscaling (min/max nodes) is a nice feature (when it works)
- Has the ability to restrict matches based on user-defined labels on the vectors
- Potentially supports incremental updates

Bad
^^^
- Index creation, even for 500 cells, took quite a long time (like 30 minutes!)
- Performance and reliability with just one node is a bit inconsistent

Ugly
^^^^
- UI is missing a lot of features
- Network configuration is needlessly complicated
- Some calls use project_id and some use project_number, some use fully-qualified paths others short names

Configuration, Creating and Deploying an Index
----------------------------------------------

Using Vertex AI Matching Engine requires several steps before matching:

#. Network Configuration (one time)
#. Create Endpoint (one time)
#. Create Index
#. Deploy Index

Constants
^^^^^^^^^

.. code-block:: bash

    # Google project, project number and region to host Vertex AI
    export PROJECT_ID="dsp-cell-annotation-service"
    export PROJECT_NUMBER=`gcloud projects describe $PROJECT_ID | grep projectNumber | cut -d"'" -f2`
    export REGION="us-central1"

    # Bucket containing CSV/AVRO of vectors to be searched
    export BUCKET_URI="gs://dsp-cell-annotation-service/demo_4m_v2/new_embeddings_for_loading/"
    export DIMENSIONS=512
    export APPROX_NEIGHBORS_COUNT=100

    # Constants, not necessary to change
    export VPC_NETWORK="ai-matching"
    export PEERING_RANGE_NAME="ann-haystack-range"
    export INDEX_ENDPOINT_NAME="casp_index_endpoint"
    export INDEX_NAME="casp_index_v1"

Network Configuration
^^^^^^^^^^^^^^^^^^^^^

.. note::

    You can skip this if network configuration has already been done. It should be done once per project.

Create the VPC Network
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    gcloud compute networks create ${VPC_NETWORK} --bgp-routing-mode=regional --subnet-mode=auto --project=${PROJECT_ID}

Add necessary firewall rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    gcloud compute firewall-rules create ${VPC_NETWORK}-allow-icmp --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow icmp

    gcloud compute firewall-rules create ${VPC_NETWORK}-allow-internal --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow all --source-ranges 10.128.0.0/9

    gcloud compute firewall-rules create ${VPC_NETWORK}-allow-rdp --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow tcp:3389

    gcloud compute firewall-rules create ${VPC_NETWORK}-allow-ssh --network ${VPC_NETWORK} --priority 65534 --project ${PROJECT_ID} --allow tcp:22

Reserve IP range
~~~~~~~~~~~~~~~~

.. code-block:: bash

    gcloud compute addresses create ${PEERING_RANGE_NAME} --global --prefix-length=16 --network=${VPC_NETWORK} --purpose=VPC_PEERING --project=${PROJECT_ID}

Set up peering with service networking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Your account must have the "Compute Network Admin" role to run the following.

.. code-block:: bash

    gcloud services vpc-peerings connect --service=servicenetworking.googleapis.com --network=${VPC_NETWORK} --ranges=${PEERING_RANGE_NAME} --project=${PROJECT_ID}

Managing Indexes
^^^^^^^^^^^^^^^^

Create Index Endpoint (to serve the index)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step takes several minutes to complete.

.. note::

    You can skip this if the endpoint has already been created. New indexes can be deployed to existing endpoints.

.. code-block:: bash

    gcloud ai index-endpoints create --display-name ${INDEX_ENDPOINT_NAME} --network projects/${PROJECT_NUMBER}/global/networks/${VPC_NETWORK} --region ${REGION} --project $PROJECT_ID

Create Index
~~~~~~~~~~~~

Creating the actual index takes a long time! (~30 minutes even for a small dataset).

.. code-block:: bash

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

This is an async operation; you will have to poll for success (the command is given by the create command above).

For example:

.. code-block:: bash

    gcloud ai operations describe 2843220864793575424 --index=7139735929568100352 --region us-central1 --project=dsp-cell-annotation-service

Deploy Index
~~~~~~~~~~~~

Deploy the index to the endpoint so it can be searched. Several non-intuitive IDs are required to run this step.

.. code-block:: bash

    # This is an identifier and a display name YOU give for this deployed index (can be the same)
    export DEPLOYED_INDEX_ID="deployed_4m_${INDEX_NAME}"
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
      --max-replica-count 2

This is an async operation; you will have to poll for success (the command is given by the create command above).

For example:

.. code-block:: bash

    gcloud ai operations describe 1574402038526115840 --index-endpoint=82032363525111808 --project $PROJECT_ID --region $REGION

Search!
-------

Searching can only be performed from compute on the same network that was configured above with the proper peering settings. The easiest way to do this is to create a Notebook instance and under the Networking configuration choose the VPC network created in the above steps (i.e., ``ai-matching`` in this example).

The DIMENSIONS, ENDPOINT_ID, and DEPLOYED_INDEX_ID variables should have the value from above.

Then from that notebook VM:

.. code-block:: python

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

Evaluating Performance
----------------------

Aspects to consider:

#. Throughput (overall matches per second)
#. Latency (response time per request)
#. Scalability (with respect to index size)
#. Accuracy
#. Cost

TBD

Cleaning Up (excluding the network setup)
-----------------------------------------

If you want to remove everything, just go in the opposite order from the above.

.. code-block:: bash

    # Undeploy Index from Endpoint
    gcloud ai index-endpoints undeploy-index ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION} --deployed-index-id=${DEPLOYED_INDEX_ID}

    # Delete Endpoint
    gcloud ai index-endpoints delete ${ENDPOINT_ID} --project ${PROJECT_ID} --region ${REGION}

    # Delete Index
    gcloud ai indexes delete ${INDEX_ID} --project ${PROJECT_ID} --region ${REGION}
