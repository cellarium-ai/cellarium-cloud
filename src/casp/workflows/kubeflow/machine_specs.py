import typing as t
from dataclasses import dataclass

from casp.workflows.kubeflow import constants

DOCKER_IMAGE_NAME_CPU = (
    "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch-pipeline:"
    "fg-embedding-job-48"
)
DOCKER_IMAGE_NAME_CUDA = (
    "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch-cuda-pipeline:"
    "fg-embedding-job-48"
)


class MachineSpecBase:
    display_name: str
    replica_count: int
    machine_type: str
    accelerator_type: t.Optional[str] = None
    accelerator_count: t.Optional[int] = None


@dataclass
class PCATrainMachineSpec(MachineSpecBase):
    display_name: str = "PCA Training"
    replica_count: int = 4
    machine_type: str = "n1-highmem-16"
    accelerator_type: t.Optional[str] = "NVIDIA_TESLA_T4"
    accelerator_count: t.Optional[int] = 4


@dataclass
class PCAEmbedMachineSpec(MachineSpecBase):
    display_name: str = "PCA Embedding"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"
    accelerator_type: str = "NVIDIA_TESLA_T4"
    accelerator_count: int = 1


@dataclass
class PCARegistryMachineSpec(MachineSpecBase):
    display_name: str = "PCA Model Registry"
    replica_count: int = 1
    machine_type: str = "n1-standard-4"


@dataclass
class PCAIndexCreateMachineSpec(MachineSpecBase):
    display_name: str = "PCA Index Create Deploy and Register"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"


@dataclass
class MeanVarStdMachineSpec(MachineSpecBase):
    display_name: str = "Mean Var Std Training"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"
    accelerator_type: str = "NVIDIA_TESLA_T4"
    accelerator_count: int = 4


@dataclass
class TDigestMachineSpec(MachineSpecBase):
    display_name: str = "TDigest Training"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"


@dataclass
class LogisticRegressionTrainMachineSpec(MachineSpecBase):
    display_name: str = "Logistic Regression Training"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"
    accelerator_type: str = "NVIDIA_TESLA_T4"
    accelerator_count: int = 1


@dataclass
class TDigestFilterFeaturesMachineSpec(MachineSpecBase):
    display_name: str = "Filter Features For TDigest"
    replica_count: int = 1
    machine_type: str = "n1-highmem-16"
    accelerator_type: str = "NVIDIA_TESLA_T4"
    accelerator_count: int = 1


component_machine_specs_map = {
    constants.PCA_TRAIN_COMPONENT_NAME: PCATrainMachineSpec,
    constants.PCA_EMBED_COMPONENT_NAME: PCAEmbedMachineSpec,
    constants.PCA_REGISTRY_COMPONENT_NAME: PCARegistryMachineSpec,
    constants.PCA_INDEX_CREATE_COMPONENT_NAME: PCAIndexCreateMachineSpec,
    constants.MEAN_VAR_STD_COMPONENT_NAME: MeanVarStdMachineSpec,
    constants.TDIGEST_COMPONENT_NAME: TDigestMachineSpec,
    constants.LOGISTIC_REGRESSION_TRAIN_COMPONENT_NAME: LogisticRegressionTrainMachineSpec,
    constants.TDIGEST_FILTER_FEATURES_COMPONENT_NAME: TDigestFilterFeaturesMachineSpec,
}
