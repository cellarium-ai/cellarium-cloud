from flask import Flask
from flask_httpauth import HTTPBasicAuth
from casp.admin.settings import settings
from casp import db

flask_app = Flask(__name__)
flask_app.config.from_object(settings)
basic_auth = HTTPBasicAuth()
db_session = db.init_db()

BASIC_AUTH_USER = {
    "username": settings.FLASK_BASIC_AUTH_USERNAME,
    "password": settings.FLASK_BASIC_AUTH_PASSWORD
}


@basic_auth.verify_password
def verify_password(username, password):
    if (
            username == BASIC_AUTH_USER.get("username")
            and password == BASIC_AUTH_USER.get("password")
    ):
        return username


with flask_app.app_context():
    db.create_all()

import casp.admin.views  # noqa
