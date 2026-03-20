.. _ci-cd:

CI/CD
=====

CAS Backend uses GitHub Actions for validation, image build, and Cloud Run deployment.

Main Workflows
--------------

- ``.github/workflows/main-workflow.yml`` runs unit-test validation on pushes
- ``.github/workflows/docker-workflow.yml`` builds and publishes service images
- ``.github/workflows/deploy-workflow.yml`` deploys service images to Cloud Run
- ``.github/workflows/integration-test-workflow.yml`` handles the integration-test flow

Validation
----------

The reusable ``.github/actions/unit-tests`` action installs Poetry, installs dependencies, runs ``make lint``, and
runs ``make test``.

Docker Build
------------

The reusable ``.github/actions/docker-build`` action selects the Dockerfile under ``deploy/docker/`` based on the
service type and pushes the resulting image to Artifact Registry.

Cloud Run Deploy
----------------

The reusable ``.github/actions/docker-deploy`` action reads service sizing from ``deploy/cloudrun/{flavor}.json`` and
deploys the selected image to Cloud Run with the configured SQL, VPC, and secret settings.

Secrets
-------

Workflow secrets and runtime environment values are documented in :doc:`secrets`.
