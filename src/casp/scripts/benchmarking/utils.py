import typing as t

from smart_open import open


def get_paths(paths: t.Union[str, t.List[str]]) -> t.List[str]:
    """
    Get a list of paths from a single path or a list of paths.

    :param paths: A single path or a list of paths.

    :return: A list of paths.
    """
    if isinstance(paths, list):
        return paths

    with open(paths, "r") as f:
        return f.read().splitlines()
