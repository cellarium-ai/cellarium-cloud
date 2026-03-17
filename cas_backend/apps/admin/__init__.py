from flask import Flask
from flask_httpauth import HTTPBasicAuth

from casp.services import db, settings

flask_app = Flask(__name__)
flask_app.config.from_object(settings)
basic_auth = HTTPBasicAuth()
db_session = db.get_db_session_maker()()


@basic_auth.verify_password
def verify_password(username: str, password: str) -> str:
    """
    Check if the provided credentials correspond to the those from the secret env variables
    :param username: Admin username
    :param password:  Admin password
    :return:
    """
    if username == settings.ADMIN_BASIC_AUTH_USER.get("username") and password == settings.ADMIN_BASIC_AUTH_USER.get(
        "password"
    ):
        return username


@flask_app.teardown_appcontext
def shutdown_session(exception=None):
    # Ensure that the session is closed at the end of each request
    db_session.close()


import casp.services.admin.views  # noqa
