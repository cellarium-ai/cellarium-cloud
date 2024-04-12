import typing as t

from casp.data_manager import BaseDataManager
from casp.services.db import models


class ModelInferenceDataManager(BaseDataManager):
    """
    Data Manager for accessing data within the scope of model inference.
    """

    def get_model_by(self, model_name: str) -> t.Optional[models.CASModel]:
        """
        Get a specific model by its unique name

        :param model_name: unique model name

        :return: CAS model object
        """
        with self.system_data_db_session_maker() as session:
            return session.query(models.CASModel).filter_by(model_name=model_name).first()
