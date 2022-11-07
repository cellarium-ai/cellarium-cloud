import torch
from casp.ml.utils import PickleMixin


class DumpManager(PickleMixin):
    def __init__(self, model, transform, running_stat):
        self.model = model
        self.transform = transform
        self.running_stat = running_stat

    @staticmethod
    def _to_device(obj, device):
        for k, attr in obj.__dict__.items():
            if isinstance(attr, torch.Tensor):
                obj.__dict__[k] = attr.to(device=device)

    def _all_attr_to_device(self, device):
        for attr_name, attr in self.__dict__.items():
            self._to_device(obj=attr, device=device)

    def to_cpu(self):
        self._all_attr_to_device(device="cpu")

    def to_cuda(self):
        self._all_attr_to_device(device="cuda")
