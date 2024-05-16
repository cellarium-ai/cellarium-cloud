import asyncio
import typing as t
from importlib import reload

from mock_alchemy.mocking import UnifiedAlchemyMagicMock
from mockito import mock, unstub, when
from sqlalchemy import orm
from sqlalchemy.orm import sessionmaker

from casp.services import db


def read_resource(resource_path: str) -> str:
    """
    Read a resource file and return its content as a string.

    Args:
        resource_path: The relative path to the resource file. (e.g. tests/unit/test_query_responses/rest_response_0.json)

    Returns:
        The content of the resource file as a string.

    """
    with open(resource_path) as file:
        return file.read()


def async_return(result: t.Any) -> asyncio.Future:
    """
    Return a result from an async function.

    Useful for mocking client calls.

    Args:
        result: The result to return.

    Returns:
        An asyncio Future that will return the result.
    """
    f = asyncio.Future()
    f.set_result(result)
    return f


def mock_sqlalchemy(data: list = []) -> None:
    """
    Mock the SQLAlchemy session for unit tests.  Be sure to call unmock_sqlalchemy() in the teardown_method of the test to reset mocks.

    See documentation for mock_alchemy at https://mock-alchemy.readthedocs.io/en/latest/user_guide/index.html#stub-data
    for more information on how to pass data in.

    Returns: The session mock object.  This can be used in the test to further configure the mock (e.g. to add objects)
    """
    mock_db_session = UnifiedAlchemyMagicMock(data=data)
    sess_maker = mock(sessionmaker)
    when(sess_maker).__call__().thenReturn(mock_db_session)

    when(orm).sessionmaker(...).thenReturn(sess_maker)

    # This is needed to make the mock work with context managers (e.g. a `with` statements)
    when(mock_db_session).__enter__(mock_db_session).thenReturn(mock_db_session)

    # Reload DB module to apply the mock
    reload(db)

    return mock_db_session


def unmock_sqlalchemy() -> None:
    """
    Unmock the SQLAlchemy session for unit tests.  This should be called in the teardown_method of the test to reset the mocks.
    """
    unstub()
    reload(db)
