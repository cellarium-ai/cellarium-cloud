import typing as t
from datetime import datetime
from uuid import UUID

from argon2 import PasswordHasher
from argon2.exceptions import Argon2Error
from sqlalchemy.exc import NoResultFound

from casp.services._auth import exceptions
from casp.services.db import get_db_session_maker, models


class OpaqueToken(t.NamedTuple):
    """
    A named tuple to represent an opaque token generated by the system.
    """

    # The user's key.  Note: this should never get stored in the database.
    key: str
    # A hashed key value
    key_hash: str


def generate_opaque_token_for_user(key_locator: UUID, key: UUID) -> OpaqueToken:
    """
    Generate an opaque token for the user. The token is a concatenation of the key locator and a key value (e.g. the secret).

    :param key_locator: A UUID key locator to find the key in the database.
    :param key: A UUID key value to authenticate the user.

    """
    ph = PasswordHasher()
    return OpaqueToken(key=f"{key_locator}.{key}", key_hash=ph.hash(str(key)))


def authenticate_user_with_opaque_token(token: str) -> models.User:
    key_locator, key = token.split(".")
    ph = PasswordHasher()
    with get_db_session_maker()() as session:
        try:
            user_key: models.UserKey = (
                session.query(models.UserKey).filter(models.UserKey.key_locator == key_locator).one()
            )
        except NoResultFound:
            raise exceptions.TokenInvalid("Token not found")
        try:
            ph.verify(user_key.key_hash, key)
        except Argon2Error:
            # Token isn't valid
            raise exceptions.TokenInvalid("Invalid Token")
        if user_key.expires < datetime.now():
            # Token is expired
            raise exceptions.TokenExpired("Token is expired")
        if not user_key.active:
            # User isn't active
            raise exceptions.TokenInactive("Token is inactive")
        if not user_key.user.active:
            # Token isn't active
            raise exceptions.TokenInactive("User is inactive")
        return user_key.user
