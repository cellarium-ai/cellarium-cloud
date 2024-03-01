Redis
=====

Summary
-------

Redis is an open source, advanced key-value store. It is often used for caching, real-time analytics, and message brokering.
`More information about Redis <http://redis.io/>`_.

Key Features:
  - Stores data in memory
  - Has all major data types like strings, hashes, lists, sets, sorted sets, etc.
  - Pub/Sub messaging

Connection
----------
To connect to redis instance, the application require the following information:
  - Host
  - Port
  - Password
  - Database

.. note::
    To connect to redis in a remote environment, the application must have access to the VPC where the redis instance is running.

Deployment
----------
In remote environment, the redis instance is deployed in a VPC as a distinct service or cluster. Locally, the redis instance can be deployed as a docker container.

Local
~~~~~
To run, you'd need to run the following commands:

.. code-block:: bash

    CONTAINER_NAME=redis-container
    PORT=6379
    PASSWORD=YourSuperSecretPassword

    docker pull redis
    docker run --name $CONTAINER_NAME -d -p $PORT:$PORT redis redis-server --requirepass $PASSWORD --port $PORT

Check if the container is running:
Installing redis-cli. On Mac OS, you can install redis (client is included) by using brew:

.. code-block:: bash

    brew install redis

Test connection:

.. code-block:: bash

    redis-cli -p $PORT -a $PASSWORD ping
    -> PONG

Remote
~~~~~~

.. note::
    The following steps are for GCP. For other cloud providers, the steps may vary.

Deploy redis instance in GCP:

.. code-block:: bash

    REDIS_INSTANCE_NAME=cellarium-cloud-redis-dev # Name of the Redis instance
    PROJECT_ID=example-project-id # GCP Project ID
    TIER=basic # basic, standard
    SIZE=4 # GB of memory
    REGION=us-central1 # Region in which the Redis instance will be created
    REDIS_VERSION=redis_7_0 # Redis version
    NETWORK=ai-matching # Network in which the Redis instance will be created
    CONNECT_MODE=PRIVATE_SERVICE_ACCESS # PRIVATE_SERVICE_ACCESS, DIRECT_PEERING
    TRANSIT_ENCRYPTION_MODE=SERVER_AUTHENTICATION # SERVER_AUTHENTICATION, DISABLED
    REDIS_VERSION=redis_7_0 # Redis version
    DISPLAY_NAME="Cellarium Cloud Redis Dev" # Display name of the Redis instance

    gcloud redis instances create $REDIS_INSTANCE_NAME \
    --project=$PROJECT_ID  \
    --tier=$TIER \
    --size=$SIZE \
    --region=$REGION \
    --redis-version=$REDIS_VERSION \
    --network=projects/dsp-cell-annotation-service/global/networks/ai-matching \
    --connect-mode=$CONNECT_MODE \
    --transit-encryption-mode=$TRANSIT_ENCRYPTION_MODE \
    --display-name=$DISPLAY_NAME \
    --enable-auth
