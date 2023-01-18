from flask import Flask
from flask_httpauth import HTTPBasicAuth

from casp import settings
from casp.services import db

flask_app = Flask(__name__)
flask_app.config.from_object(settings)
basic_auth = HTTPBasicAuth()
db_session = db.init_db()


@basic_auth.verify_password
def verify_password(username, password):
    if username == settings.ADMIN_BASIC_AUTH_USER.get("username") and password == settings.ADMIN_BASIC_AUTH_USER.get(
        "password"
    ):
        return username


with flask_app.app_context():
    pass

import casp.services.admin.views  # noqa
