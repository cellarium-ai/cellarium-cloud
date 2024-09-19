import datetime
import typing as t

from packaging.version import Version

from casp.services import constants
from casp.services.api import schemas
from casp.services.api.data_manager import CellariumGeneralDataManager
from casp.services.api.services import exceptions
from casp.services.db import models


class CellariumGeneralService:
    """
    Cellarium General Service. Service for managing data and access to general Cellarium Cloud information.
    """

    def __init__(self, cellarium_general_dm: CellariumGeneralDataManager = None):
        self.cellarium_general_dm = cellarium_general_dm or CellariumGeneralDataManager()

    def get_application_info(self) -> schemas.ApplicationInfo:
        """
        Get application information

        :return: Object with CAS application information
        """
        return self.cellarium_general_dm.get_application_info()

    def validate_client_version(self, client_version: str) -> bool:
        """
        Check whether the client version is new enough to work with the server

        :param client_version: Client version string
        """
        try:
            return Version(client_version) >= Version(constants.MIN_CLIENT_VERSION)
        except ValueError:
            raise exceptions.InvalidClientVersionException(client_version=client_version)

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
            weekly_quota=user.cell_quota,
            remaining_weekly_quota=self.cellarium_general_dm.get_remaining_quota_for_user(user=user),
            quota_reset_date=quota_reset_date,
            lifetime_quota=user.lifetime_cell_quota,
            remaining_lifetime_quota=(
                None if user.lifetime_cell_quota is None else user.lifetime_cell_quota - user.total_cells_processed
            ),
        )
