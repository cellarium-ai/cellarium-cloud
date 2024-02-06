.. _ci-cd:

CI/CD
=====
Currently repository git actions are configured to run tests and update docker images on every push.
To change this behavior, edit `.github/workflows/main.yml` file. Main workflow is separated into actions (for test and docker build).

Test
****
Test job is run through tox. Tox is configured to run tests with pytest and python code formatting and linting tools (black, isort, flake8).

Docker Update
*************
All docker images are stored in `casp/services/deploy/<docker_image_name>`. Action job creates a matrix of docker images for different purposes.

Docker Images
^^^^^^^^^^^^^
- `casp/services/deploy/Dockerfile.pytorch` - PyTorch image with installed Cellarium Cloud services. Used to deploy services that require PyTorch without GPU support.
- `casp/services/deploy/Dockerfile.pytorch_cuda` - PyTorch image with installed Cellarium Cloud services. Used to deploy services that require PyTorch with GPU support.

Secrets
^^^^^^^
Docker update actions requires GitHub secrets to be set. Please refer to :doc:`secrets` for secret details.
