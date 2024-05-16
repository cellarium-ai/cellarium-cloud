import json
import typing as t
from base64 import b64decode, b64encode
from datetime import datetime, timedelta
from unittest import mock
from uuid import UUID, uuid4

import jwt
import pytest
from fastapi.security import HTTPAuthorizationCredentials
from mockito import unstub, when
from parameterized import parameterized
from starlette_context import context, request_cycle_context

from casp.services import _auth, settings
from casp.services._auth import OpaqueToken, exceptions, jwt_token, opaque_token
from casp.services.api.dependencies.auth import TokenTypes, authenticate_user
from casp.services.constants import ContextKeys
from casp.services.db import models
from tests.unit.test_utils import mock_sqlalchemy, unmock_sqlalchemy

TEST_HASHING_KEY: str = "test"


class TestAuth:
    """
    Test various auth mechanisms in the application.
    """

    def setup_method(self) -> None:
        self.mock_db_session = mock_sqlalchemy()

    def teardown_method(self) -> None:
        unstub()
        unmock_sqlalchemy()

    @parameterized.expand(
        [
            # Encoded JWT token with payload {"a": "b"} and hashing key "a"
            ("eyJhbGciOiJIUzI1NiJ9.eyJhIjoiYiJ9.bdAkj_C_GHDuLuFTBxOGTtT8mxPXLFvrIlHGX9sgq7o", TokenTypes.JWT, None),
            ("00000000-0000-0000-0000-000000000000.11111111-1111-1111-1111-111111111111", TokenTypes.OPAQUE, None),
            (
                "Bearer 00000000-0000-0000-0000-000000000000.11111111-1111-1111-1111-111111111111",
                TokenTypes.OPAQUE,
                None,
            ),
            ("foo", TokenTypes.OPAQUE, exceptions.TokenInvalid),
            ("a.b.c.d", TokenTypes.OPAQUE, exceptions.TokenInvalid),
            ("", TokenTypes.OPAQUE, exceptions.TokenInvalid),
        ]
    )
    def test_token_type_detection(
        self, token: str, expected_type: TokenTypes, expected_exception: t.Optional[Exception]
    ) -> None:
        if expected_exception:
            with pytest.raises(expected_exception):
                TokenTypes.detect_type(token)
        else:
            assert TokenTypes.detect_type(token) == expected_type

    @parameterized.expand(
        [
            (
                HTTPAuthorizationCredentials(
                    credentials=jwt_token.generate_jwt_for_user(user_id=1, token_ttl=100, hashing_key=TEST_HASHING_KEY),
                    scheme="bearer",
                ),
                None,
            ),
            (
                HTTPAuthorizationCredentials(
                    credentials=opaque_token.generate_opaque_token_for_user(key_locator=uuid4(), key=uuid4()).key,
                    scheme="bearer",
                ),
                None,
            ),
            (HTTPAuthorizationCredentials(credentials="foo", scheme="bearer"), exceptions.TokenInvalid),
            (HTTPAuthorizationCredentials(credentials="a.b.c.d", scheme="bearer"), exceptions.TokenInvalid),
            (None, exceptions.TokenInvalid),
        ]
    )
    @pytest.mark.asyncio
    async def test_authenticate_user(
        self, auth_token_scheme: HTTPAuthorizationCredentials, expected_exception: t.Optional[Exception]
    ):
        """
        Test the top level authentication function.
        """

        user_id: int = 1
        user = models.User(id=user_id)

        if expected_exception:
            with pytest.raises(expected_exception):
                await authenticate_user(auth_token_scheme)
        else:
            when(_auth).authenticate_user_with_opaque_token(...).thenReturn(user)
            when(_auth).authenticate_user_with_jwt(...).thenReturn(user)

            with request_cycle_context() as _:
                assert await authenticate_user(auth_token_scheme) == user
                assert context[ContextKeys.user] == user

    def test_generate_jwt_for_user(self) -> None:
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)
        parsed_token = jwt.decode(token, TEST_HASHING_KEY, algorithms=[settings.JWT_HASHING_ALGORITHM])

        assert parsed_token["user_id"] == user_id
        assert datetime.fromisoformat(parsed_token["expiration"]) > datetime.now()

    def test_authenticate_user_jwt(self) -> None:
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)

        # Since we are just doing a primary key lookup we can just mock the object (as opposed to the query)
        user = models.User(id=user_id, active=True)
        self.mock_db_session.add(user)

        assert jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY) == user

    def test_authenticate_user_jwt_not_found(self) -> None:
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)

        with pytest.raises(exceptions.TokenInvalid):
            jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY)

    def test_authenticate_user_jwt_not_active(self) -> None:
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)

        # Since we are just doing a primary key lookup we can just mock the object (as opposed to the query)
        user = models.User(id=user_id, active=False)
        self.mock_db_session.add(user)

        with pytest.raises(exceptions.TokenInactive):
            jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY)

    def test_authenticate_user_jwt_expired(self) -> None:
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=-10, hashing_key=TEST_HASHING_KEY)

        # Since we are just doing a primary key lookup we can just mock the object (as opposed to the query)
        user = models.User(id=user_id, active=False)
        self.mock_db_session.add(user)

        with pytest.raises(exceptions.TokenExpired):
            jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY)

    def test_authenticate_user_with_jwt_with_bad_hash_key(self) -> None:
        """
        Test that bad actors can't authenticate with a key generated with a different hashing key.
        """
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)
        with pytest.raises(exceptions.TokenInvalid):
            jwt_token.authenticate_user_with_jwt(token, hashing_key="bad_key")

    def test_authenticate_user_with_jwt_with_tampered_key(self) -> None:
        """
        Test that bad actors can't authenticate with a key generated that has been tampered with
        """
        user_id: int = 1
        token = jwt_token.generate_jwt_for_user(user_id=user_id, token_ttl=100, hashing_key=TEST_HASHING_KEY)

        # replace the payload of the JWT with a different user_id
        header, payload, signature = token.split(".")
        parsed_payload = json.loads(b64decode(payload + "=="))  # extra padding needed for b64decode
        parsed_payload["user_id"] = 2
        tampered_payload = b64encode(json.dumps(parsed_payload).encode())
        tampered_jwt = f"{header}.{tampered_payload}.{signature}"

        with pytest.raises(exceptions.TokenInvalid):
            jwt_token.authenticate_user_with_jwt(tampered_jwt, hashing_key=TEST_HASHING_KEY)

    def test_authenticate_user_with_jwt_with_(self) -> None:
        """
        Test that bad actors can't authenticate with a key generated with a different hashing key.
        """
        user_id: int = 1
        payload = {"user_id": user_id, "foo": "bar"}
        token = jwt.encode(payload=payload, key=TEST_HASHING_KEY, algorithm=settings.JWT_HASHING_ALGORITHM)
        with pytest.raises(exceptions.TokenInvalid):
            jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY)

    def test_authenticate_user_with_jwt_with_(self) -> None:
        """
        Test that bad actors can't authenticate with a key generated with a different hashing key.
        """
        user_id: int = 1
        payload = {"user_id": user_id, "foo": "bar"}
        token = jwt.encode(payload=payload, key=TEST_HASHING_KEY, algorithm=settings.JWT_HASHING_ALGORITHM)
        with pytest.raises(exceptions.TokenInvalid):
            jwt_token.authenticate_user_with_jwt(token, hashing_key=TEST_HASHING_KEY)

    @parameterized.expand(
        [
            # Success case
            (
                models.User(id=10, active=True),
                20,
                uuid4(),
                uuid4(),
                False,
                True,
                datetime.now() + timedelta(hours=1),
                False,
                None,
            ),
            # Inactive user
            (
                models.User(id=11, active=False),
                21,
                uuid4(),
                uuid4(),
                False,
                True,
                datetime.now() + timedelta(hours=1),
                False,
                exceptions.TokenInactive,
            ),
            # Inactive user key
            (
                models.User(id=12, active=True),
                22,
                uuid4(),
                uuid4(),
                False,
                False,
                datetime.now() + timedelta(hours=1),
                False,
                exceptions.TokenInactive,
            ),
            # Invalid user key
            (
                models.User(id=13, active=True),
                23,
                uuid4(),
                uuid4(),
                True,
                True,
                datetime.now() + timedelta(hours=1),
                False,
                exceptions.TokenInvalid,
            ),
            # Expired token
            (
                models.User(id=14, active=True),
                24,
                uuid4(),
                uuid4(),
                False,
                True,
                datetime.now() - timedelta(hours=1),
                False,
                exceptions.TokenExpired,
            ),
            # User not found
            (
                models.User(id=15, active=True),
                25,
                uuid4(),
                uuid4(),
                False,
                True,
                datetime.now() + timedelta(hours=1),
                True,
                exceptions.TokenInvalid,
            ),
        ]
    )
    def test_authenticate_user_with_opaque_token(
        self,
        user: models.User,
        key_id: int,
        key_locator: UUID,
        key: UUID,
        make_bad_key: bool,
        key_active: bool,
        expires: datetime,
        bad_key_locator: bool,
        expected_exception: Exception,
    ) -> None:
        token: OpaqueToken = opaque_token.generate_opaque_token_for_user(key_locator=key_locator, key=key)

        key_hash = token.key_hash
        user_key = models.UserKey(
            id=key_id,
            key_locator=str(key_locator),
            key_hash=key_hash,
            user_id=user.id,
            user=user,
            expires=expires,
            active=key_active,
        )

        # Since we want to mock filtering we need to mock the query as opposed to the objects
        self.mock_db_session = mock_sqlalchemy(
            data=[
                (
                    [
                        mock.call.query(models.UserKey),
                        mock.call.filter(
                            models.UserKey.key_locator == str(key_locator if not bad_key_locator else uuid4())
                        ),
                    ],
                    [user_key],
                ),
            ]
        )

        key_to_use = token.key
        if make_bad_key:
            bad_token: OpaqueToken = opaque_token.generate_opaque_token_for_user(key_locator=key_locator, key=uuid4())
            key_to_use = bad_token.key
        if expected_exception is None:
            assert opaque_token.authenticate_user_with_opaque_token(key_to_use) == user
        else:
            with pytest.raises(expected_exception):
                opaque_token.authenticate_user_with_opaque_token(key_to_use)
