Deploying in Cloud
==================

CAS Backend services are deployed on `Google Cloud Run <https://cloud.google.com/run/docs>`_ using the shared GitHub
Actions workflows and the deployment assets under ``deploy/``.

Deployment Assets
-----------------

- ``deploy/docker/`` contains service Dockerfiles
- ``deploy/cloudrun/`` contains flavor-specific Cloud Run sizing
- ``deploy/scripts/`` contains supporting deployment scripts
- ``.github/actions/docker-build`` and ``.github/actions/docker-deploy`` contain the reusable CI/CD actions

Service Images
--------------

Current deployable service images:

- compute: ``deploy/docker/Dockerfile.compute``
- admin: ``deploy/docker/Dockerfile.admin``

Deployment Flow
---------------

1. Export runtime dependencies into ``deploy/requirements.txt.lock`` with ``make requirements``.
2. Build and push the service image.
3. Deploy the image to Cloud Run using the configured flavor, service account, VPC connector, SQL instance, and secret.

The repository GitHub workflows automate this flow for standard deployments.

Configuration
-------------

Cloud Run sizing and concurrency settings are selected from the JSON files in ``deploy/cloudrun/``.
Secrets are mounted into the container at ``/app/settings/.env``.

Related Documentation
---------------------

- :doc:`readme_modules/cloudrun`
- :doc:`cicd`
- :doc:`secrets`

.. toctree::
   :maxdepth: 1

   readme_modules/cloudrun
