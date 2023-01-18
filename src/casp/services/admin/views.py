import tempfile

from flask import Response, redirect, request, send_file
from flask_admin import Admin, AdminIndexView, expose
from flask_admin.contrib.sqla import ModelView
from flask_admin.model.template import EndpointLinkRowAction
from werkzeug.exceptions import HTTPException

from casp.services import _auth
from casp.services.admin import basic_auth, flask_app
from casp.services.db import db_session, models


class AuthException(HTTPException):
    def __init__(self, message):
        super().__init__(
            message,
            Response(
                "You could not be authenticated. Please refresh the page.",
                401,
                {"WWW-Authenticate": 'Basic realm="Login Required"'},
            ),
        )


class BasicHTTPAuthMixin:
    @staticmethod
    def is_accessible():
        auth_obj = basic_auth.get_auth()
        password = basic_auth.get_auth_password(auth_obj)

        if not basic_auth.authenticate(auth_obj, password):
            raise AuthException("Not authenticated.")

        return True

    @staticmethod
    def inaccessible_callback(name, **kwargs):
        return redirect(basic_auth.challenge())


class CASAdminIndexView(BasicHTTPAuthMixin, AdminIndexView):
    pass


class CASAdminModelView(BasicHTTPAuthMixin, ModelView):
    pass


class UserAdminView(CASAdminModelView):
    column_extra_row_actions = [
        EndpointLinkRowAction("glyphicon glyphicon-asterisk", ".generate_secret_key"),
    ]
    column_list = ("email", "is_active", "cas_request_count")
    form_widget_args = {"cas_request_count": {"disabled": True}, "cas_scRNA_cells_processed": {"disabled": True}}

    _token = None
    _page_refreshed = False

    @staticmethod
    def _create_token_file(token) -> tempfile.NamedTemporaryFile:
        temp = tempfile.NamedTemporaryFile()
        temp.name = "api_token.txt"
        temp.write(token.encode())
        temp.seek(0)
        return temp

    @expose("/generate-secret-key", methods=("GET",))
    def generate_secret_key(self):
        token = _auth.generate_jwt_for_user(int(request.args["id"]))
        temp = self._create_token_file(token)
        return send_file(temp, as_attachment=True, download_name=temp.name)


admin = Admin(
    flask_app,
    name="CAS Admin",
    template_mode="bootstrap3",
    index_view=CASAdminIndexView(url="/", template="admin/main_page.html"),
)
admin.add_view(UserAdminView(models.User, db_session, name="User"))
