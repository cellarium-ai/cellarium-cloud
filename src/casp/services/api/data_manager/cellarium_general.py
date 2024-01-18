import typing as t

from casp.data_manager import BaseDataManager, sql
from casp.services import settings
from casp.services.api import schemas
from casp.services.api.data_manager import exceptions
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

    @staticmethod
    def get_models_all() -> t.List[models.CASModel]:
        """
        Retrieve query with all available models

        :return: List of CAS models
        """
        return models.CASModel.query.all()

    @staticmethod
    def get_models_non_admin() -> t.List[models.CASModel]:
        """
        Retrieve query with available models only for non-admin use (admin_use_only=False)

        :return: List of CAS models
        """
        return models.CASModel.query.filter_by(admin_use_only=False).all()

    @staticmethod
    def get_model_by_name(model_name: str) -> models.CASModel:
        """
        Retrieve CAS model by its system name

        :param model_name: ML Model name

        :raises: NotFound if model is not found

        :return: CAS ML model object from the database
        """
        model = models.CASModel.query.filter_by(model_name=model_name).first()

        if model is None:
            raise exceptions.NotFound(f"Model {model_name} not found in the database")

        return model

    @classmethod
    def get_index_for_model(cls, model_name: str) -> models.CASMatchingEngineIndex:
        """
        Retrieve CAS model index by its system name

        :param model_name: ML Model name

        :return: CAS ML model index object from the database
        """
        return cls.get_model_by_name(model_name=model_name).cas_matching_engine

    def increment_user_cells_processed(self, user: models.User, number_of_cells: int) -> None:
        """
        Increment `requests_processed` and `cells_processed user stats`.

        :param user: A SQLAlchemy user model
        :param number_of_cells: number of cells to increment by
        """
        user.requests_processed += 1
        user.cells_processed += number_of_cells
        self.system_data_db_session.commit()
