import tempfile

from flask import Response, request, send_file, redirect
from flask_admin import Admin, AdminIndexView, expose
from flask_admin.contrib.sqla import ModelView
from flask_admin.model.template import EndpointLinkRowAction
from werkzeug.exceptions import HTTPException

from casp.services import _auth, utils
from casp.services.admin import basic_auth, flask_app
from casp.services.db import db_session, models, ops


class AuthException(HTTPException):
    def __init__(self, message):
        super().__init__(
            message,
            Response(
                response="You could not be authenticated. Please refresh the page.",
                status=401,
                headers={"WWW-Authenticate": 'Basic realm="Login Required"'},
            ),
        )


class BasicHTTPAuthMixin:
    """
    Basic HTTP Auth Mixin
    Implements HTTP Basic Auth in 2 methods that are used by `Flask-Admin` views:
    `is_accessible`, `inaccessible_callback`
    """

    @staticmethod
    def is_accessible() -> bool:
        """
        Authorize the user request by HTTP Basic Auth rule. Throw an AuthException if the request wasn't
        authorized
        :return: True if the request was authorized
        """
        auth_obj = basic_auth.get_auth()
        password = basic_auth.get_auth_password(auth_obj)

        if not basic_auth.authenticate(auth_obj, password):
            raise AuthException("Not authenticated.")

        return True


class CellariumCloudAdminIndexView(BasicHTTPAuthMixin, AdminIndexView):
    """
    Index Page for the project Admin dashboard.
    Inherits a Basic HTTP protection rule from `BasicHTTPAuthMixin`
    """


class CellariumCloudAdminModelView(BasicHTTPAuthMixin, ModelView):
    """
    Model View class for Admin Dashboard.
    Inherits a Basic HTTP protection rule from `BasicHTTPAuthMixin`
    """

    pass


class UserAdminView(CellariumCloudAdminModelView):
    """
    User Admin View. Has a custom method `generate-secret-key` which is used to generate a
    JWT token to let the user call the authorization protected methods.
    """

    column_extra_row_actions = [
        EndpointLinkRowAction("glyphicon glyphicon-asterisk", ".generate_secret_key"),
    ]
    column_list = ("email", "is_admin", "is_active", "requests_processed", "cells_processed")
    column_editable_list = ("is_admin",)
    form_widget_args = {"requests_processed": {"disabled": True}, "cells_processed": {"disabled": True}}

    @staticmethod
    def _create_token_file(token) -> tempfile.NamedTemporaryFile:
        """
        Generate a temporary file with JWT token
        :param token: JWT Token
        :return: Temporary file with a JWT token string
        """
        temp = tempfile.NamedTemporaryFile()
        temp.name = "api_token.txt"
        temp.write(token.encode())
        temp.seek(0)
        return temp

    @expose("/generate-secret-key", methods=("GET",))
    def generate_secret_key(self) -> Response:
        """
        Exposes an API for Flask Admin for generating secret key (JWT token) and
        downloading it as a .txt file
        :return: A HTTP response that makes a web browser downloading the file
        """
        token = _auth.generate_jwt_for_user(int(request.args["id"]))
        temp = self._create_token_file(token)
        return send_file(temp, as_attachment=True, download_name=temp.name)


class CASModelAdminView(CellariumCloudAdminModelView):
    column_extra_row_actions = [
        EndpointLinkRowAction("glyphicon glyphicon-chevron-up", ".set_default_model"),
    ]
    column_list = (
        "system_name",
        "model_file_path",
        "embedding_dimension",
        "admin_use_only",
        "schema_name",
        "is_default_model",
        "bq_cell_info_table_fqn",
        "bq_temp_table_dataset",
        "created_date",
    )
    column_descriptions = {
        "system_name": (
            "A system name that is used, must be unique, lowercase. "
            "No spaces, must end with a character or number. \nExample: cas-pca-001."
        ),
        "model_file_path": "Filepath in the GCS storage bucket with the dumped model.",
        "embedding_dimension": "Model embedding output dimension.",
        "admin_use_only": (
            "Flag switching the access to this model to all the users. "
            "If false, only admin users can access the model endpoint. "
            "Set this to false when model is tested and well benchmarked."
        ),
        "schema_name": "Schema name that was used in data to build the model.",
        "is_default_model": (
            "Flag showing to CAS client which model schema to use as default. Only one model could be a default model."
        ),
        "bq_cell_info_table_fqn": (
            "Cell info table that was used to create extract the model was trained on. "
            "This would be used for joining cell info for annotations"
        ),
        "bq_temp_table_dataset": (
            "Dataset where CAS would allocate temporary tables for storing annotation request data."
        ),
        "created_date": "Datetime when this record has been created. Differs from when model was trained.",
    }
    column_editable_list = ("admin_use_only",)
    form_columns = (
        "system_name",
        "model_file_path",
        "embedding_dimension",
        "schema_name",
        "is_default_model",
        "admin_use_only",
        "bq_cell_info_table_fqn",
        "bq_temp_table_dataset",
        "created_date",
    )
    form_widget_args = {"created_date": {"disabled": True}, "is_default_model": {"disabled": True}}

    @expose("/set-default-model", methods=("GET",))
    def set_default_model(self) -> Response:
        """
        Exposes an API for Flask Admin for generating secret key (JWT token) and
        downloading it as a .txt file
        :return: A HTTP response that makes a web browser downloading the file
        """
        model_id = int(request.args["id"])
        ops.set_default_model_by(model_id=model_id)
        return redirect("/casmodel/")


class CASMatchingEngineAdminView(CellariumCloudAdminModelView):
    column_list = ("system_name", "embedding_dimension", "endpoint_id", "deployed_index_id", "model")
    column_descriptions = {
        "system_name": (
            "A system name that is used, must be unique, lowercase. "
            "No spaces, must end with a character or number. \nExample: cas-pca-001-matching-engine-index."
        ),
        "endpoint_id": "Endpoint ID that is used in GCP in Vertex AI",
        "deployed_index_id": "Deployed Index ID that is used in GCP in Vertex AI",
    }


admin = Admin(
    flask_app,
    name="Cellarium Cloud Admin",
    template_mode="bootstrap3",
    index_view=CellariumCloudAdminIndexView(url="/", template="admin/main_page.html"),
)
admin.add_view(UserAdminView(models.User, db_session, name="User"))
admin.add_view(CASModelAdminView(models.CASModel, db_session, name="CASModel"))
admin.add_view(CASMatchingEngineAdminView(models.CASMatchingEngineIndex, db_session, name="MatchingEngine"))
