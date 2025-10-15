import asyncio

from casp.services.api import schemas
from casp.services.api.services.annotation_engines.model_interpretation_engine import (
    ModelInterpretationStrategyInterface,
)
from casp.services.model_inference.schemas import ModelInferenceOutputBase


class ModelInterpretationEngine:
    """
    Model interpretation engine responsible for summarizing model output.
    """

    def __init__(self, strategy: ModelInterpretationStrategyInterface):
        self.strategy = strategy

    def summarize(self, model_output: ModelInferenceOutputBase) -> schemas.QueryAnnotationAbstractType:
        """
        Summarize query neighbor context using the specified strategy.

        :param model_output: Output of the model inference

        :return: A list of post-processed model responses.
        """
        return self.strategy.summarize(model_output=model_output)

    async def summarize_async(self, model_output: ModelInferenceOutputBase) -> schemas.QueryAnnotationAbstractType:
        """
        Call :meth:`summarize` in a separate thread to avoid blocking the event loop.

        :param model_output: Output of the model inference

        :return: A list of post-processed model responses.
        """
        return await asyncio.to_thread(self.summarize, model_output)
