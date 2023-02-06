import logging
from datetime import datetime, timedelta

import jwt
from sqlalchemy.exc import IntegrityError

from casp.services import settings
from casp.services._auth import exceptions
from casp.services.db import models


def generate_jwt_for_user(user_id: int, token_ttl: int = settings.JWT_DEFAULT_TOKEN_TTL) -> str:
    """
    Create a JWT token for the user. JWT contains 2 parameters in a payload: `user_id`
    and `expiration`
    :param user_id: A user id to reference them in a db
    :param token_ttl: Token time to live
    :return: An JWT token for the further authentication
    """
    payload = {"user_id": user_id, "expiration": str(datetime.now() + timedelta(seconds=token_ttl))}
    return jwt.encode(payload=payload, key=settings.SECURITY_PASSWORD_SALT, algorithm=settings.JWT_HASHING_ALGORITHM)


def authenticate_user_with_jwt(token: str) -> models.User:
    """
    Authenticate a user with JWT. Throw `TokenInvalid` if token couldn't be parsed
    or `TokenExpired` if it is parsable but expired.
    :param token: JWT token
    :return: A SQLAlchemy User model from the db
    """
    try:
        payload = jwt.decode(
            jwt=token, key=settings.SECURITY_PASSWORD_SALT, algorithms=[settings.JWT_HASHING_ALGORITHM]
        )
    except jwt.PyJWTError as e:
        logging.error(e)
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
    except IntegrityError as e:
        logging.error(e)
        raise exceptions.TokenInvalid
