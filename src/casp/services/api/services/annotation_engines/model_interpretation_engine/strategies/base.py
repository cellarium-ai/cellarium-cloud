import typing as t

from casp.services.api import schemas
from casp.services.model_inference.schemas import ModelInferenceOutputBase


class ModelInterpretationStrategyInterface:
    def summarize(self, model_output: ModelInferenceOutputBase) -> t.List[schemas.QueryCellAnnotationAbstract]:
        """
        Define the method signature for the consensus strategy, which is responsible for summarizing query neighbor

        :param model_output: Output from model inference service

        :return: A list of summarized annotations per querying cell.
        """
        raise NotImplementedError
