
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
