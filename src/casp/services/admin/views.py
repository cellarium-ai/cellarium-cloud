import logging
import tempfile
from uuid import UUID, uuid4

import typing_extensions as tx
from flask import Response, flash, redirect, request, send_file
from flask_admin import Admin, AdminIndexView, expose
from flask_admin.babel import gettext
from flask_admin.contrib.sqla import ModelView
from flask_admin.contrib.sqla.fields import InlineModelFormList
from flask_admin.contrib.sqla.form import InlineModelConverter
from flask_admin.form import BaseForm, FormOpts, RenderTemplateWidget
from flask_admin.helpers import get_redirect_target
from flask_admin.model.helpers import get_mdict_item_or_list
from flask_admin.model.template import EndpointLinkRowAction, LinkRowAction
from werkzeug.exceptions import HTTPException

from casp.services import _auth, settings
from casp.services.admin import basic_auth, db_session, flask_app
from casp.services.db import models, ops
from casp.services.utils.email_utils import EmailSender

email_sender = EmailSender(sendgrid_key=settings.SENDGRID_API_KEY, from_address=settings.FROM_ADDRESS)

logger = logging.getLogger(__name__)


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


class UserKeyInlineFieldListWidget(RenderTemplateWidget):
    """
    Widget that uses a custom template for the user key inline field list
    """

    def __init__(self):
        super().__init__("admin/inline_list_user_key.html")


class UserKeylineModelFormList(InlineModelFormList):
    """
    Used in the chain of classes needed to use a custom template
    """

    widget = UserKeyInlineFieldListWidget()


class UserKeyInlineModelConverter(InlineModelConverter):
    """
    Used in the chain of classes needed to use a custom template
    """

    inline_field_list_type = UserKeylineModelFormList


class UserAdminView(CellariumCloudAdminModelView):
    """
    User Admin View. Has a custom method `generate-secret-key` which is used to generate a
    JWT token to let the user call the authorization protected methods.
    """

    column_extra_row_actions = [
        EndpointLinkRowAction("glyphicon glyphicon-asterisk", ".generate_secret_key"),
        LinkRowAction("glyphicon glyphicon-duplicate", "clone?id={row_id}"),
    ]
    column_list = (
        "email",
        "is_admin",
        "active",
        "cell_quota",
        "lifetime_cell_quota",
        "total_requests_processed",
        "total_cells_processed",
        "created_at",
    )
    column_editable_list = ("is_admin",)
    form_excluded_columns = ("requests_processed", "cells_processed")

    inline_models = [
        (
            models.UserKey,
            dict(form_label="User Keys", form_excluded_columns=["created_date", "key_hash", "key_locator"]),
        )
    ]
    inline_model_form_converter = UserKeyInlineModelConverter

    # Use as a mechanism to pass keys from the on_model_change to the after_model_change callbacks since we
    # do not want to store the actual user key in the database
    key_cache: dict[str, str] = {}

    @tx.override
    def on_model_change(self, _: BaseForm, model: models.User, is_created: bool):

        if is_created and len(model.user_keys) == 0:
            # Expiration and id are created when inserting into the database
            key_locator: UUID = uuid4()
            key: UUID = uuid4()
            token: _auth.OpaqueToken = _auth.generate_opaque_token_for_user(key_locator, key)

            model.user_keys.append(
                models.UserKey(
                    key_locator=str(key_locator),
                    key_hash=token.key_hash,
                    active=True,
                    user=model,
                )
            )

            # Cache the key to be shown during the after_model_change invocation
            self.key_cache[str(key_locator)] = token.key
        else:
            # Find which key(s) were created and cache them to be shown during the after_model_change invocation
            for user_key in model.user_keys:
                if user_key.key_locator is None:
                    key_locator: UUID = uuid4()
                    key: UUID = uuid4()
                    token: _auth.OpaqueToken = _auth.generate_opaque_token_for_user(key_locator, key)

                    # Set the keylocator and keyhash to store in the database
                    user_key.key_locator = str(key_locator)
                    user_key.key_hash = token.key_hash

                    self.key_cache[str(key_locator)] = token.key

    @tx.override
    def after_model_change(self, form: BaseForm, model: models.User, is_created: bool):
        if is_created and len(model.user_keys) == 1:
            # Notifies the user that the key was created.
            # Note: the brackets are used for the javascript to parse the key out of the message
            try:
                key = self.key_cache[model.user_keys[0].key_locator]
                # Renove the key from the cache
                del self.key_cache[model.user_keys[0].key_locator]
                email_sender.send_welcome_email(email=model.email, key=key)
            except Exception:
                logger.error("Error notifying user", exc_info=True)
                flash(gettext(f"Error notifying user {model.email}. Please contact them manually."), "error")
        else:
            # Find which key(s) were created and print the key to the console as a one shot deal
            for user_key in model.user_keys:
                if user_key.key_locator in self.key_cache:
                    try:
                        key = self.key_cache[user_key.key_locator]
                        del self.key_cache[user_key.key_locator]
                        email_sender.send_new_key_email(email=model.email, key=key)
                    except Exception:
                        logger.error("Error notifying user", exc_info=True)
                        flash(gettext(f"Error notifying user {model.email}. Please contact them manually."), "error")
        return super().after_model_change(form, model, is_created)

    @staticmethod
    def _create_token_file(token: str) -> tempfile.NamedTemporaryFile:
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
        "description",
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
        "description": "A more verbose description of the model, since the name is not always self-explanatory.",
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
        "description",
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
    column_list = (
        "index_name",
        "description",
        "embedding_dimension",
        "endpoint_id",
        "deployed_index_id",
        "num_neighbors",
        "model",
    )
    column_descriptions = {
        "index_name": (
            "A name that is used to identify the index, must be unique, lowercase. "
            "No spaces, must end with a character or number. \nExample: cas-pca-001-matching-engine-index."
        ),
        "description": "A more verbose description of the index, since the name is not always self-explanatory.",
        "endpoint_id": "Endpoint ID that is used in GCP in Vertex AI",
        "deployed_index_id": "Deployed Index ID that is used in GCP in Vertex AI",
        "num_neighbors": "Number of neighbors that is used for an approximate neighbors search",
    }
    column_extra_row_actions = [
        LinkRowAction("glyphicon glyphicon-duplicate", "clone?id={row_id}"),
    ]


def shorten_value_formatter(view, context, model, name) -> str:
    """
    Shorten the value to 20 characters plus ellipsis or return empty string if None
    """
    value = getattr(model, name)

    if value is None:
        return ""

    return value[:20] + "..." if len(value) > 20 else value


class CellInfoAdminView(CellariumCloudAdminModelView):
    column_list = (
        "cas_cell_index",
        "original_cell_id",
        "cell_type",
        "cell_type_ontology_term_id",
        "assay",
        "development_stage",
        "tissue",
        "disease",
        "cas_ingest.dataset_id",
    )
    column_labels = {"cas_ingest.dataset_id": "Dataset ID"}
    column_filters = (
        "cas_cell_index",
        "cas_ingest_id",
        "original_cell_id",
        "cell_type",
        "cell_type_ontology_term_id",
        "assay",
        "development_stage",
        "tissue",
        "disease",
        "cas_ingest.dataset_id",
    )
    column_select_related_list = ("cas_ingest",)
    column_formatters = {"original_cell_id": shorten_value_formatter}
    # Do not allow to edit or delete the records, neither create new ones. This is a read-only view.
    can_edit = False
    can_delete = False
    can_create = False

    page_size = 50


class CellFeatureInfoAdminView(CellariumCloudAdminModelView):
    column_list = (
        "cas_feature_index",
        "original_feature_id",
        "feature_name",
        "feature_biotype",
        "feature_is_filtered",
        "feature_reference",
    )
    column_filters = (
        "cas_feature_index",
        "original_feature_id",
        "feature_name",
        "feature_biotype",
        "feature_is_filtered",
        "feature_reference",
    )
    column_formatters = {"original_feature_id": shorten_value_formatter}
    page_size = 50


class CellIngestInfoAdminView(CellariumCloudAdminModelView):
    column_list = ("cas_ingest_id", "dataset_id", "ingest_timestamp")
    column_filters = ("cas_ingest_id", "dataset_id", "ingest_timestamp")


admin = Admin(
    flask_app,
    name="Cellarium Cloud Admin",
    template_mode="bootstrap3",
    index_view=CellariumCloudAdminIndexView(url="/", template="admin/main_page.html"),
)
admin.add_view(UserAdminView(models.User, db_session, name="User", category="Users"))
admin.add_view(CASModelAdminView(models.CASModel, db_session, name="CASModel", category="ML Management"))
admin.add_view(
    CASMatchingEngineAdminView(
        models.CASMatchingEngineIndex, db_session, name="MatchingEngine", category="ML Management"
    )
)
admin.add_view(CellInfoAdminView(models.CellInfo, db_session, name="CellInfo", category="Cell Data Management"))
admin.add_view(
    CellFeatureInfoAdminView(models.FeatureInfo, db_session, name="CellFeature", category="Cell Data Management")
)
admin.add_view(
    CellIngestInfoAdminView(models.CellIngestInfo, db_session, name="CellIngestInfo", category="Cell Data Management")
)
