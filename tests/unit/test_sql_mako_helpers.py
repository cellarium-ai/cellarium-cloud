"""
Test cases for the functions in the 'mako_helpers' module of the 'casp.data_manager' package.

These test cases cover the following functions:
- _string_value_processor
- _bool_value_processor
- select_clause
- where_clause

"""

import pytest

from casp.data_manager.sql import mako_helpers


def test_string_value_processor() -> None:
    """
    Test the _string_value_processor function.

    This function tests the processing of string values by adding single quotes around them.

    Test cases:
    - It should add single quotes around the provided string.
    - It should return an empty string if an empty string is provided.
    """
    assert mako_helpers._string_value_processor("test") == "'test'"
    assert mako_helpers._string_value_processor("hello world") == "'hello world'"
    assert mako_helpers._string_value_processor("") == "''"


def test_bool_value_processor() -> None:
    """
    Test the _bool_value_processor function.

    This function tests the processing of boolean values by converting them to uppercase strings.

    Test cases:
    - It should convert True to 'TRUE'.
    - It should convert False to 'FALSE'.
    """
    assert mako_helpers._bool_value_processor(True) == "TRUE"
    assert mako_helpers._bool_value_processor(False) == "FALSE"


# Test cases for select_clause
def test_select_clause() -> None:
    """
    Test the select_clause function.

    This function tests the construction of a SQL SELECT clause from a list of columns.

    Test cases:
    - When an empty list is provided, it should select all columns ('select *').
    - When a list of columns is provided, it should construct the SELECT clause correctly.
    """
    assert mako_helpers.select([]) == "*"
    assert mako_helpers.select(["column1", "column2"]) == "select column1, column2"
    assert mako_helpers.select(["column"]) == "select column"


def test_where_clause() -> None:
    """
    Test the where_clause function.

    This function tests the construction of a SQL WHERE clause from filter criteria.

    Test cases:
    - When no filters or an empty dictionary is provided, it should return an empty string.
    - When 'equals', 'in', 'not equals', 'not in' filter types are used, it should construct the WHERE clause correctly.
    - When an unsupported filter type is provided, it should raise a ValueError.
    - When string and boolean values are used in 'equals' filter, it should construct the WHERE clause correctly.
    """
    # Test with no filters
    assert mako_helpers.where({}) == ""

    # Test 'equals' filter
    assert mako_helpers.where({"organism__eq": "Homo sapiens"}) == "where\n    organism = 'Homo sapiens'"

    # Test 'in' filter
    assert mako_helpers.where({"cell_type__in": ["T cell", "neuron"]}) == "where\n    cell_type in ('T cell', 'neuron')"

    # Test 'not equals' filter
    assert mako_helpers.where({"organism__not_eq": "Homo sapiens"}) == "where\n    organism != 'Homo sapiens'"

    # Test 'not in' filter
    assert (
        mako_helpers.where({"cell_type__not_in": ["T cell", "neuron"]})
        == "where\n    cell_type not in ('T cell', 'neuron')"
    )

    # Test unsupported filter type
    with pytest.raises(ValueError):
        mako_helpers.where({"column__unsupported": "value"})

    # Test string and boolean values in 'equals' filter
    assert mako_helpers.where({"name__eq": "John"}) == "where\n    name = 'John'"
    assert mako_helpers.where({"is_active__eq": True}) == "where\n    is_active = TRUE"
