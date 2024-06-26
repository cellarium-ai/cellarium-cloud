import datetime
import typing as t

from casp.services.api import schemas
from casp.services.api.data_manager import CellariumGeneralDataManager
from casp.services.api.services import exceptions
from casp.services.db import models


class CellariumGeneralService:
    """
    Cellarium General Service. Service for managing data and access to general Cellarium Cloud information.
    """

    def __init__(self):
        self.cellarium_general_dm = CellariumGeneralDataManager()

    def get_application_info(self) -> schemas.ApplicationInfo:
        """
        Get application information

        :return: Object with CAS application information
        """
        return self.cellarium_general_dm.get_application_info()

    def feedback_opt_out(self, user: models.User) -> models.User:
        """
        Opt out user from feedback requests

        :param user: User object to opt out
        """
        return self.cellarium_general_dm.feedback_opt_out(user=user)

    def get_feature_schemas(self) -> t.List[schemas.FeatureSchemaInfo]:
        """
        Get all available feature schemas

        :return: List of gene schema objects
        """
        return self.cellarium_general_dm.get_feature_schemas()

    def get_feature_schema_by(self, schema_name: str) -> t.List[str]:
        """
        Get a specific feature schema by its unique name

        :param schema_name: unique feature schema name

        :return: List of features in a correct order
        """
        return self.cellarium_general_dm.get_feature_schema_by(schema_name=schema_name)

    def get_model_list_for_user(self, user: models.User) -> t.List[models.CASModel]:
        """
        Get query with available models for a specific user based on their permissions.

        :param user: User object to check permissions for

        :return: List of CAS models
        """
        if user.is_admin:
            return self.cellarium_general_dm.get_models_all()

        return self.cellarium_general_dm.get_models_non_admin()

    def get_model_for_user(self, model_name: str, user: models.User) -> models.CASModel:
        """
        Get CAS model by its system name for a specific user based on their permissions.

        :param model_name: Name of the model
        :param user: User object to check permissions for

        :return: CAS ML model object from the database
        """
        model = self.cellarium_general_dm.get_model_by_name(model_name=model_name)

        if model.admin_use_only and not user.is_admin:
            raise exceptions.AccessDeniedError(
                f"{model_name} model is not available. Please reach out to the Cellarium team for more information."
            )
        return model

    def get_quota_for_user(self, user: models.User) -> schemas.UserQuota:
        """
        Get user quota information

        :param user: User object to check quota for

        :return: User quota information
        """

        days_to_add = datetime.datetime.now().weekday() - user.created_at.weekday()
        if days_to_add <= 0:
            days_to_add += 7
        quota_reset_date = datetime.datetime.today().replace(
            hour=0, minute=0, second=0, microsecond=0
        ) + datetime.timedelta(days=days_to_add)

        return schemas.UserQuota(
            user_id=user.id,
            quota=user.cell_quota,
            remaining_quota=self.cellarium_general_dm.get_remaining_quota_for_user(user=user),
            quota_reset_date=quota_reset_date,
        )
