import datetime
import typing as t

import sqlalchemy as sa
from starlette_context import context

from casp.data_manager import BaseDataManager, sql
from casp.services import settings
from casp.services.api import schemas
from casp.services.api.data_manager import exceptions
from casp.services.constants import ContextKeys
from casp.services.db import models


class CellariumGeneralDataManager(BaseDataManager):
    """
    Data Access Object for managing data within the scope of general Cellarium Cloud information.
    """

    # Directories for SQL templates
    SQL_TEMPLATE_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_templates/cellarium_general"

    # SQL template file paths
    SQL_GET_ALL_GENE_SCHEMAS = f"{SQL_TEMPLATE_DIR}/get_all_gene_schemas.sql.mako"
    SQL_GET_SCHEMA_BY_NAME = f"{SQL_TEMPLATE_DIR}/get_schema_by_name.sql.mako"

    # SQL query directory
    SQL_QUERY_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_queries"

    # SQL queries from files
    with open(f"{SQL_QUERY_DIR}/get_cells_processed_this_week_for_user.sql", "r") as f:
        SQL_GET_CELLS_PROCESSED_THIS_WEEK = f.read()



    @staticmethod
    def get_application_info() -> schemas.ApplicationInfo:
        """
        :return: Object with CAS application information
        """
        return schemas.ApplicationInfo(
            default_feature_schema=settings.DEFAULT_FEATURE_SCHEMA, application_version=settings.APP_VERSION
        )

    def get_feature_schemas(self) -> t.List[schemas.FeatureSchemaInfo]:
        """
        :return: List of gene schema objects
        """
        sql_template_data = sql.TemplateData(project=self.project)
        sql_query = sql.render(template_path=self.SQL_GET_ALL_GENE_SCHEMAS, template_data=sql_template_data)

        query_result = self.block_coo_matrix_db_client.query(sql_query).result()
        return [schemas.FeatureSchemaInfo(schema_name=row["table_name"]) for row in query_result]

    def get_feature_schema_by(self, schema_name: str) -> t.List[str]:
        """
        Get a specific feature schema by its unique name and return a list of features in a correct order

        :param schema_name: unique feature schema name

        :return: List of features in a correct order
        """
        sql_template_data = sql.TemplateData(project=self.project, schema_name=schema_name)
        sql_query = sql.render(template_path=self.SQL_GET_SCHEMA_BY_NAME, template_data=sql_template_data)
        query_result = self.block_coo_matrix_db_client.query(sql_query).result()
        return [row["feature_name"] for row in query_result]

    def get_models_all(self) -> t.List[models.CASModel]:
        """
        Retrieve query with all available models

        :return: List of CAS models
        """
        with self.system_data_db_session_maker() as session:
            return session.query(models.CASModel).all()

    def get_models_non_admin(self) -> t.List[models.CASModel]:
        """
        Retrieve query with available models only for non-admin use (admin_use_only=False)

        :return: List of CAS models
        """
        with self.system_data_db_session_maker() as session:
            return session.query(models.CASModel).filter_by(admin_use_only=False).all()

    def get_model_by_name(self, model_name: str) -> models.CASModel:
        """
        Retrieve CAS model by its system name

        :param model_name: ML Model name

        :raises: NotFound if model is not found

        :return: CAS ML model object from the database
        """
        with self.system_data_db_session_maker() as session:
            model = session.query(models.CASModel).filter_by(model_name=model_name).first()

            if model is None:
                raise exceptions.NotFound(f"Model {model_name} not found in the database")

            return model

    def get_index_for_model(self, model_name: str) -> models.CASMatchingEngineIndex:
        """
        Retrieve CAS model index by its system name

        :param model_name: ML Model name

        :return: CAS ML model index object from the database
        """

        with self.system_data_db_session_maker() as session:
            model = session.query(models.CASModel).filter_by(model_name=model_name).first()

            if model is None:
                raise exceptions.NotFound(f"Model {model_name} not found in the database")

            return model.cas_matching_engine

    def log_user_activity(
        self, user_id: int, model_name: str, method: str, cell_count: int, event: models.UserActivityEvent
    ) -> None:
        """
        Log the information to the DB about a request from a user

        :param user_id: The id for a user
        :param model_name: Model name that user is using
        :param method: Method that user is using
        :param cell_count: Number of cells processed in this request
        """
        user_activity = models.UserActivity(
            user_id=user_id,
            request_id=context[ContextKeys.sentry_trace_id],
            model_name=model_name,
            method=method,
            cell_count=cell_count,
            event=event,
        )

        with self.system_data_db_session_maker.begin() as session:
            session.add(user_activity)

    def get_remaining_quota_for_user(self, user: models.User) -> int:
        """
        Get remaining cell processing quota for a specific user (i.e. their quota minus cells processed this week)

        :param user: User object to check quota for

        :return: Remaining cell processing quota
        """
        # Return the remaining quota, or 0 if it's negative (which might happen if the user exceeds
        # their quota, either because we allow them to, or because of any weirdness in the
        # simultaneous processing of cells)
        return max(user.cell_quota - self.get_cells_processed_this_week_for_user(user=user), 0)

    def get_cells_processed_this_week_for_user(self, user: models.User) -> int:
        """
        Get number of cells processed by a specific user this week

        :param user: User object to check cells processed for

        :return: Number of cells processed this week
        """
        # Start of week for a user is based on the day of the week they signed up, so we need to
        # use the day of the week they signed up to calculate the start of the week
        current_weekday = datetime.datetime.today().weekday()
        start_of_the_week_for_user = user.created_at.weekday()
        days_to_subtract = current_weekday - start_of_the_week_for_user
        if days_to_subtract < 0:
            days_to_subtract += 7
        start_of_week = datetime.datetime.today().replace(
            hour=0, minute=0, second=0, microsecond=0
        ) - datetime.timedelta(days=days_to_subtract)

        with self.system_data_db_session_maker() as session:
            cell_count_sum_query = sa.sql.text(self.SQL_GET_CELLS_PROCESSED_THIS_WEEK)
            return session.execute(cell_count_sum_query, {"id": user.id, "start_of_week": start_of_week}).first()[0]
