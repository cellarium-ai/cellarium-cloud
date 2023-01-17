from datetime import datetime, timedelta
import jwt
from sqlalchemy.exc import IntegrityError
from casp.admin import settings
from casp.db import models
from casp.auth import exceptions


def generate_token_for_user(user_id: int) -> str:
    payload = {"user_id": user_id, "expiration": str(datetime.now() + timedelta(days=60))}
    return jwt.encode(payload, settings.SECURITY_PASSWORD_SALT, algorithm="HS256")


def authenticate_token(token: str) -> models.User:
    try:
        payload = jwt.decode(token, settings.SECURITY_PASSWORD_SALT, algorithms=["HS256"])
    except jwt.PyJWTError:
        raise exceptions.TokenInvalid

    if "user_id" not in payload or "expiration" not in payload:
        raise exceptions.TokenInvalid

    expiration_date_str = payload["expiration"]
    user_id = payload["user_id"]
    expiration_date = datetime.fromisoformat(expiration_date_str)

    if datetime.now() >= expiration_date:
        raise exceptions.TokenExpired
    try:
        return models.User.query.get(user_id)
    except IntegrityError:
        raise exceptions.TokenInvalid
