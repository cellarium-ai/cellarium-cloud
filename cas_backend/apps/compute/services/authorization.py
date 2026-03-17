import typing as t

from casp.services.api.data_manager import CellariumGeneralDataManager
from casp.services.api.data_manager import exceptions as dm_exc
from casp.services.api.services import exceptions
from casp.services.db import models


class Authorizer:
    """
    Service for authorizing users to access certain resources.
    """

    def __init__(self, cellarium_general_dm: t.Optional[CellariumGeneralDataManager] = None):
        self.cellarium_general_dm = cellarium_general_dm or CellariumGeneralDataManager()

    def authorize_model_for_user(self, user: models.User, model_name: str) -> models.CASModel:
        """
        Verify the user's authorization to use a specific model. If the user is not authorized, raise an exception.

        :param user: User object to check permissions for.
        :param model_name: Model name to check permissions for.

        :return: Model object if the model exists and if the user is authorized.
        """
        try:
            model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)
        except dm_exc.NotFound as e:
            raise exceptions.InvalidInputError(str(e))
        if model.admin_use_only and not user.is_admin:
            raise exceptions.AccessDeniedError(
                f"{model_name} model is not available. Please reach out to the Cellarium team for more information."
            )
        return model
