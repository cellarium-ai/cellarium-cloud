import tempfile
from werkzeug.exceptions import HTTPException
from flask import request, send_file, Response, redirect
from flask_admin.contrib.sqla import ModelView
from flask_admin import Admin, AdminIndexView
from flask_admin import expose
from flask_admin.model.template import EndpointLinkRowAction

from casp.db import models, db_session
from casp import auth
from casp.admin import flask_app, basic_auth


class AuthException(HTTPException):
    def __init__(self, message):
        super().__init__(message, Response(
            "You could not be authenticated. Please refresh the page.", 401,
            {'WWW-Authenticate': 'Basic realm="Login Required"'}))


class BasicHTTPAuthMixin:
    @staticmethod
    def is_accessible():
        auth = basic_auth.get_auth()
        password = basic_auth.get_auth_password(auth)

        if not basic_auth.authenticate(auth, password):
            raise AuthException('Not authenticated.')
        else:
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
        EndpointLinkRowAction("glyphicon glyphicon-copy", ".generate_secret_key"),
    ]
    column_list = ("email", "is_active",)

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
        token = auth.generate_token_for_user(int(request.args["id"]))
        temp = self._create_token_file(token)
        return send_file(temp, as_attachment=True, download_name=temp.name)


class APITokenAdminView(CASAdminModelView):
    column_formatters = {"key": lambda v, c, m, p: f"{m.key[:3]}***{m.key[-3:]}"}


admin = Admin(
    flask_app,
    name="Cell Annotation Service Admin",
    template_mode="bootstrap3",
    index_view=CASAdminIndexView(url="/"),
)
admin.add_view(UserAdminView(models.User, db_session, name="User"))
