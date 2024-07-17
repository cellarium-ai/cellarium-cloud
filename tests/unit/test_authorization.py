import re

import pytest
from mockito import mock, unstub, when

from casp.services.api.data_manager import exceptions as dm_exc
from casp.services.api.services import Authorizer, exceptions
from casp.services.db import models

USER_ADMIN = models.User(id=1, is_admin=True)
USER_NON_ADMIN = models.User(id=2, is_admin=False)

MODEL = models.CASModel(id=1, model_name="model_name", admin_use_only=False)
MODEL_ADMIN_ONLY = models.CASModel(id=2, model_name="admin_only_model", admin_use_only=True)
NON_EXISTENT_MODEL = models.CASModel(id=3, model_name="non_existent_model", admin_use_only=True)


class TestAuthorization:
    def setup_method(self) -> None:
        self.authorizer = Authorizer(cellarium_general_dm=mock())

    def teardown_method(self) -> None:
        unstub()

    def test_non_admin_user_can_accss_model(self) -> None:
        self.__mock_model_queries(model=MODEL)
        assert self.authorizer.authorize_model_for_user(user=USER_NON_ADMIN, model_name=MODEL.model_name) == MODEL

    def test_admin_user_can_accss_admin_model(self) -> None:
        self.__mock_model_queries(model=MODEL_ADMIN_ONLY)
        assert (
            self.authorizer.authorize_model_for_user(user=USER_ADMIN, model_name=MODEL_ADMIN_ONLY.model_name)
            == MODEL_ADMIN_ONLY
        )

    def test_non_admin_can_not_access_admin_model(self) -> None:
        self.__mock_model_queries(model=MODEL_ADMIN_ONLY)
        with pytest.raises(exceptions.AccessDeniedError, match=re.escape("admin_only_model model is not available.")):
            self.authorizer.authorize_model_for_user(user=USER_NON_ADMIN, model_name=MODEL_ADMIN_ONLY.model_name)

    def test_model_does_not_exist(self) -> None:
        when(self.authorizer.cellarium_general_dm).get_model_by_name(model_name="non_existent_model").thenRaise(
            dm_exc.NotFound()
        )
        with pytest.raises(exceptions.InvalidInputError):
            self.authorizer.authorize_model_for_user(user=USER_NON_ADMIN, model_name=NON_EXISTENT_MODEL.model_name)

    def __mock_model_queries(
        self,
        model: models.CASModel,
    ) -> None:
        """
        Mock call postgres that fetched the model

        :param model: The model to mock.
        """

        when(self.authorizer.cellarium_general_dm).get_model_by_name(model_name=model.model_name).thenReturn(model)
