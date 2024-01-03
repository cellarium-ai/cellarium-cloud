TRAIN_DISPLAY_NAME = "PCA Training"
TRAIN_REPLICA_COUNT = 1
TRAIN_MACHINE_TYPE = "n1-highmem-16"
TRAIN_ACCELERATOR_TYPE = "NVIDIA_TESLA_T4"
TRAIN_ACCELERATOR_COUNT = 4

EMBED_DISPLAY_NAME = "PCA Embedding"
EMBED_REPLICA_COUNT = 1
EMBED_MACHINE_TYPE = "n1-highmem-16"
EMBED_ACCELERATOR_TYPE = "NVIDIA_TESLA_T4"
EMBED_ACCELERATOR_COUNT = 4

REGISTRY_DISPLAY_NAME = "PCA Model Registry"
REGISTRY_REPLICA_COUNT = 1
REGISTRY_MACHINE_TYPE = "n1-standard-4"
# REGISTRY_NETWORK_NAME = "ai-matching"
REGISTRY_NETWORK_NAME = "projects/350868384795/global/networks/ai-matching"

INDEX_CREATE_DISPLAY_NAME = "PCA Index Create Deploy and Register"

DOCKER_IMAGE_NAME_CPU = (
    "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:fg-embedding-job-23"
)
DOCKER_IMAGE_NAME_CUDA = (
    "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch-cuda:fg-embedding-job-23"
)
