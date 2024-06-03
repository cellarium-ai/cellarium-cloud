Deployment Flavors
==================

This directory should contain various deployment configurations for the sizes of services used by the deployment GitHub action.
These are JSON files named `[flavor].json` and the contents should look like:

.. code-block:: json

    {
      "cas-admin": {
        "cpu": "1000m",
        "memory": "512Mi",
        "max-instances": 100,
        "min-instances": 0,
        "concurrency": 80

      },
      "cas-model": {
        "cpu": "4000m",
        "memory": "16Gi",
        "max-instances": 200,
        "min-instances": 0,
        "concurrency": 10
      },
      "cas-api": {
        "cpu": "1000m",
        "memory": "4Gi",
        "max-instances": 200,
        "min-instances": 0,
        "concurrency": 20
      }
    }

With the appropriate values for the service.