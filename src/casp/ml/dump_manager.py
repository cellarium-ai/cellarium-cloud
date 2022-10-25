from casp.ml.utils import PickleMixin


class DumpManager(PickleMixin):
    def __init__(self, model, transform, running_stat):
        self.model = model
        self.transform = transform
        self.running_stat = running_stat
