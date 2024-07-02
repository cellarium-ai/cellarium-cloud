.. _ci-cd:

CI/CD
=====
Currently repository GitHub actions are configured to run unit tests and update docker images on every push.
To change this behavior, edit `.github/workflows/main.yml` file. Main workflow is separated into actions (for test and docker build).

Test
****

Unit Test
^^^^^^^^^
Unit tests are run through tox. Tox is configured to run tests with pytest and also run python code formatting and linting tools (black, isort, flake8).

Integration Test
^^^^^^^^^^^^^^^^
Integration tests are run using tox, but are a separate workflow (`.github/workflows/integration-test-workflow.yml`) from the unit tests.
This workflow has to be run manually using either the GitHub web UI or the GitHub CLI tool.  The integration tests are defined in the cellarium-cas
repository.  Running the integration test workflow deploys the services for the specified version of cellarium-cloud to Cloud Run and runs the
integration tests defined in the specified version of cellarium-cas against the deployed services.

Docker
******

Docker Update
^^^^^^^^^^^^^
All docker images are stored in `casp/services/deploy/<docker_image_name>`. Action job creates a matrix of docker images for different purposes.

Docker Deploy
^^^^^^^^^^^^^
Docker deploy action is used to deploy the docker images to Cloud Run. The action is defined in `.github/workflows/deploy-workflow.yml`.

Docker Images
^^^^^^^^^^^^^
- `casp/services/deploy/Dockerfile.pytorch` - PyTorch image with installed Cellarium Cloud services. Used to deploy services that require PyTorch without GPU support.
- `casp/services/deploy/Dockerfile.pytorch_cuda` - PyTorch image with installed Cellarium Cloud services. Used to deploy services that require PyTorch with GPU support.

Secrets
*******
Docker update and integration test actions require GitHub secrets to be set. Please refer to :doc:`secrets` for secret details.
