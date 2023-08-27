import typing as t

import pytest

from casp.cell_data_manager import sql
from casp.cell_data_manager.sql.validation import exceptions, template_data_validator


@pytest.fixture
def valid_template_data() -> "sql.TemplateData":
    return sql.TemplateData(
        project="test-project",
        dataset="test_dataset",
        select=["test_column", "test_column_2", "test_column_3"],
        filters={"test_column_3__eq": 5, "test_column__in": ["foo", "buz", "ex"]},
    )


@pytest.fixture
def valid_filter_statements() -> t.Dict:
    return {"test_column_3__eq": 5, "test_column__in": ["foo", "buz", "ex"]}


@pytest.fixture
def invalid_filter_statements() -> t.Dict:
    return {"test_column__llk": "some_value"}


def test_validate_sql_filter_valid():
    """
    Test the validation of valid SQL filters.

    This test checks the validation of SQL filters with valid filter types (e.g., "__eq", "__in").
    It ensures that the validation function does not raise any exceptions.

    """
    template_data_validator.validate_sql_filter("column_name__eq", 42)
    template_data_validator.validate_sql_filter("column_name__in", [1, 2, 3])


def test_validate_sql_filter_invalid():
    """
    Test the validation of invalid SQL filters.

    This test checks the validation of SQL filters with an invalid filter type.
    It ensures that the validation function raises a `TemplateDataValidationError` exception.

    """
    with pytest.raises(exceptions.TemplateDataValidationError):
        template_data_validator.validate_sql_filter("column_name__invalid", 42)


def test_validate_filter_statements(valid_filter_statements):
    """
    Test the validation of valid filter statements.

    This test validates a list of valid filter statements using the `validate_filters` function.
    It ensures that the validation function does not raise any exceptions.

    """
    template_data_validator.validate_filters(valid_filter_statements)


def test_validate_filter_statements_invalid(invalid_filter_statements):
    """
    Test the validation of invalid filter statements.

    This test validates a list of invalid filter statements using the `validate_filters` function.
    It ensures that the validation function raises a `TemplateDataValidationError` exception.

    """
    with pytest.raises(exceptions.TemplateDataValidationError):
        template_data_validator.validate_filters(invalid_filter_statements)


def test_render_valid_template(valid_template_data):
    """
    Test rendering a valid SQL template.

    This test renders a valid SQL template using the `sql.render` function.
    It compares the rendered SQL query with the expected SQL query from a file.

    """
    template_path = "tests/unit/test_sql_templates/valid_sql_template.sql.mako"
    expected_sql_path = "tests/unit/test_sql_templates/valid_sql_template_expected.sql"

    # Turn off SQL Query validator because it requires a connection to BigQuery
    rendered_sql = sql.render(template_path, valid_template_data, sql_query_validator_on=False)

    with open(expected_sql_path) as file:
        expected_sql = file.read()

    assert rendered_sql == expected_sql, "Rendered SQL query does not correspond to the expected one"
