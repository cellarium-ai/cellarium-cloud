import tempfile

from flask import Response, flash, redirect, request, send_file
from flask_admin import Admin, AdminIndexView, expose
from flask_admin.babel import gettext
from flask_admin.contrib.sqla import ModelView
from flask_admin.form import FormOpts
from flask_admin.helpers import get_redirect_target
from flask_admin.model.helpers import get_mdict_item_or_list
from flask_admin.model.template import EndpointLinkRowAction, LinkRowAction
from werkzeug.exceptions import HTTPException

from casp.services import _auth
from casp.services.admin import basic_auth, flask_app
from casp.services.db import get_db_session_maker, models, ops

db_session = get_db_session_maker()()


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

    @expose("/clone/", methods=("GET", "POST"))
    def clone_view(self):
        """
        Clone model view.
        Derived from flask_admin.model.BaseModelView.create_view and flask_admin.model.BaseModelView.edit_view
        """
        return_url = get_redirect_target() or self.get_url(".index_view")

        if not self.can_create:
            return redirect(return_url)

        id = get_mdict_item_or_list(request.args, "id")
        if id is None:
            return redirect(return_url)

        model = self.get_one(id)

        if model is None:
            flash(gettext("Record does not exist."), "error")
            return redirect(return_url)

        form = self.create_form(obj=model)
        if not hasattr(form, "_validated_ruleset") or not form._validated_ruleset:
            self._validate_form_instance(ruleset=self._form_edit_rules, form=form)

        if self.validate_form(form):
            if self.create_model(form):
                flash(gettext("Record was successfully created."), "success")
                if "_add_another" in request.form:
                    return redirect(self.get_url(".create_view", url=return_url))
                elif "_continue_editing" in request.form:
                    return redirect(self.get_url(".edit_view", id=self.get_pk_value(model)))
                else:
                    # save button
                    return redirect(self.get_save_return_url(model, is_created=False))

        if request.method == "GET" or form.errors:
            self.on_form_prefill(form, id)

        form_opts = FormOpts(widget_args=self.form_widget_args, form_rules=self._form_edit_rules)

        if self.create_modal and request.args.get("modal"):
            template = self.create_modal_template
        else:
            template = self.create_template

        return self.render(template, model=model, form=form, form_opts=form_opts, return_url=return_url)


class UserAdminView(CellariumCloudAdminModelView):
    """
    User Admin View. Has a custom method `generate-secret-key` which is used to generate a
    JWT token to let the user call the authorization protected methods.
    """

    column_extra_row_actions = [
        EndpointLinkRowAction("glyphicon glyphicon-asterisk", ".generate_secret_key"),
        LinkRowAction("glyphicon glyphicon-duplicate", "clone?id={row_id}"),
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
        LinkRowAction("glyphicon glyphicon-duplicate", "clone?id={row_id}"),
    ]
    column_list = (
        "model_name",
        "model_file_path",
        "embedding_dimension",
        "admin_use_only",
        "schema_name",
        "is_default_model",
        "bq_dataset_name",
        "created_date",
    )
    column_descriptions = {
        "model_name": (
            "A name that is used, must be unique, lowercase. "
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
        "bq_dataset_name": (
            "Bigquery dataset name that is used to store cell information which were used to train the model. "
        ),
        "created_date": "Datetime when this record has been created. Differs from when model was trained.",
    }
    column_editable_list = ("admin_use_only",)
    form_columns = (
        "model_name",
        "model_file_path",
        "embedding_dimension",
        "schema_name",
        "is_default_model",
        "admin_use_only",
        "bq_dataset_name",
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
    column_list = ("index_name", "embedding_dimension", "endpoint_id", "deployed_index_id", "num_neighbors", "model")
    column_descriptions = {
        "index_name": (
            "A name that is used to identify the index, must be unique, lowercase. "
            "No spaces, must end with a character or number. \nExample: cas-pca-001-matching-engine-index."
        ),
        "endpoint_id": "Endpoint ID that is used in GCP in Vertex AI",
        "deployed_index_id": "Deployed Index ID that is used in GCP in Vertex AI",
        "num_neighbors": "Number of neighbors that is used for an approximate neighbors search",
    }
    column_extra_row_actions = [
        LinkRowAction("glyphicon glyphicon-duplicate", "clone?id={row_id}"),
    ]


admin = Admin(
    flask_app,
    name="Cellarium Cloud Admin",
    template_mode="bootstrap3",
    index_view=CellariumCloudAdminIndexView(url="/", template="admin/main_page.html"),
)
admin.add_view(UserAdminView(models.User, db_session, name="User"))
admin.add_view(CASModelAdminView(models.CASModel, db_session, name="CASModel"))
admin.add_view(CASMatchingEngineAdminView(models.CASMatchingEngineIndex, db_session, name="MatchingEngine"))
