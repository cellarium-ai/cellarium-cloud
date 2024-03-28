import asyncio
import typing as t


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
