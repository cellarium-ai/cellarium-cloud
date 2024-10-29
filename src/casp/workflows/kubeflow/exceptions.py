class KubeflowException(Exception):
    pass


class KubeflowComponentException(KubeflowException):
    pass


class MachineSpecsDoesntExist(KubeflowException):
    pass
