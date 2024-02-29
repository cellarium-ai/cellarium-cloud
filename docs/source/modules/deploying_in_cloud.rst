Deploying in Cloud
==================
Cellarium Cloud services are deployed on `Google Cloud Run <https://cloud.google.com/run/docs>`_.
These services are containerized and built and deployed through GitHub Actions. :ref:`More info <ci-cd>`.

.. code-block::

    gcloud run deploy <service_name> --project $PROJECT_ID --image $IMAGE_NAME

Instructions on how to update a particular service are provided in the respective service's documentation. :ref:`More info <services>`.

Once the services are deployed, they will appear in the `Google Cloud Run Console <https://console.cloud.google.com/run>`_
as a new instance, or an existing one will be updated if the service name already exists.

Accessing the Services
++++++++++++++++++++++

To get access to the deployed services, you can check the URL in the
`Google Cloud Run Console <https://console.cloud.google.com/run>`_. The dashboard provides logs, metrics, and error
reporting.