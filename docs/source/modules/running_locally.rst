Running Locally
===============
To run Cellarium Cloud services locally, you will need to have the following installed:

- Python 3.10
- `Local Postgres Database cluster <https://www.docker.com/blog/how-to-use-the-postgres-docker-official-image>`_
- Project Environment variables in a `src/settings/.env` file. :ref:`More info <Project Secrets>`.
- Project dependencies installed. ``pip install -r requirements.txt``
- `src` directory added to your ``PYTHONPATH`` environment variable. E.g. ``export PYTHONPATH=$PYTHONPATH:/path/to/cellarium-cloud/src``


Once you have the above installed, you can run the services locally using the following command:

.. code-block:: bash

    python src/casp/services/<service_name>/main.py

To check the API methods that exist and their documentation, you can visit `API docs page <http://localhost:8000/api/docs>`_

Most of the methods will reuqire you to be authenticated. To do this, you'd need to deploy Admin service:

.. code-block:: bash

    python src/casp/services/admin/server.py


Once the Admin service is running, you can visit `Admin Dashboard <http://127.0.0.1:5000>`_ to create a user and token.

To run commands that require accessing the vector search API using gRPC, you will need to create a tunnel.  To do this you
must:
- log into the Broad non-split VPN
- log into the gcloud cli with the command ``gcloud auth login``

You can then run the following command to create a tunnel:

.. code-block:: bash

    gcloud compute ssh --zone "us-central1-a" "bastion" --project "dsp-cell-annotation-service" -- -NL 10000:localhost:10000

In local development mode, the user agent header sent to the vector search API contains information on what vector API endpoint
to route to.

The bastion server is a proxy server that routes requests to the vector search API and is configured using nginx with a reverse
proxy that is configured in file ``/etc/nginx/sites-available/vertex``:

.. code-block:: nginx

    map $http_user_agent $forward_ip {
    "~target_ip\:([^\s]+)" $1;
    default "no-ip-to-forward-to";
    }

    server {
        listen 10000 http2;

        location / {
            grpc_pass grpc://$forward_ip:10000;
            include proxy_params;
        }
    }
