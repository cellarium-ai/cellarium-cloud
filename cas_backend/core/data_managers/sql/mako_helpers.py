import typing as t

from casp.data_manager.sql import exceptions
from casp.data_manager.sql.constants import CAS_CELL_INFO_REQUIRED_COLUMN_NAMES, ComparisonOperators
from casp.data_manager.sql.validation.template_data_validator import validate_column_name


def _string_value_processor(v: str) -> str:
    """
    Add single quotes around the provided string, preparing it for SQL syntax.

    :param v: The string to be processed.
    :return: The processed string with single quotes added.
    """
    if "'" in v:
        return f'"{v}"'

    return f"'{v}'"


def _numeric_value_processor(v: int | float) -> str:
    """
    Convert input value into string as is, preparing it for SQL syntax.

    :param v: The value to be processed.
    :return: The processed value in string format.
    """
    return str(v)


def _bool_value_processor(v: bool) -> str:
    """
    Convert the provided boolean to a string and capitalize it, making it suitable for SQL syntax.

    :param v: The boolean value to be processed.
    :return: The processed boolean in string format and uppercase."""
    return str(v).upper()


def _list_value_processor(v: list) -> str:
    """
    Convert the provided list of values to a SQL compatible string format based on the type of the first element in the
    list.

    :param v: The list of values to be processed.

    :return: The processed list of values in string format.
    """
    if isinstance(v[0], str):
        return f"({', '.join([_string_value_processor(x) for x in v])})"
    elif isinstance(v[0], int):
        return f"({', '.join([_numeric_value_processor(x) for x in v])})"


def _process_equality_operator_filter_value(filter_value: t.Any, filter_type: str) -> str:
    """
    Process the provided filter value for an equality operator.

    :param filter_value: The value to be processed.

    :return: The processed value in string format.
    """
    # TODO: Add support for date and time values
    if isinstance(filter_value, str):
        return _string_value_processor(v=filter_value)
    elif isinstance(filter_value, bool):
        return _bool_value_processor(v=filter_value)
    elif isinstance(filter_value, (int, float)):
        return _numeric_value_processor(v=filter_value)
    else:
        raise exceptions.UnsupportedSQLTypeException(
            f"Unsupported type {type(filter_value)} for `{filter_type}`. " f"Value should be a primitive python type."
        )


def _process_relational_operator_filter_value(filter_value: t.Any, filter_type: str) -> str:
    """
    Process the provided filter value for a relational operator.

    :param filter_value: The value to be processed.

    :return: The processed value in string format.
    """
    if isinstance(filter_value, (int, float)):
        return _numeric_value_processor(v=filter_value)
    else:
        raise exceptions.UnsupportedSQLTypeException(
            f"Unsupported type {type(filter_value)} of value for filter {filter_type}. "
            f"`{filter_type}` filter type supports only integers and floats."
        )


def _process_membership_operator_filter_value(filter_value: t.Any, filter_type: str) -> str:
    """
    Process the provided filter value for a membership operator.

    :param filter_value: The value to be processed.

    :return: The processed value in string format.
    """
    if isinstance(filter_value, (list, tuple)):
        return _list_value_processor(v=filter_value)
    else:
        raise exceptions.UnsupportedSQLTypeException(
            f"Unsupported type {type(filter_type)} for `{filter_type}`. "
            f"Value should be a list of numbers or strings."
        )


FILTER_VALUE_PROCESS_HANDLERS_MAPPING = {
    ComparisonOperators.EQUAL: _process_equality_operator_filter_value,
    ComparisonOperators.NOT_EQUAL: _process_equality_operator_filter_value,
    ComparisonOperators.IN: _process_membership_operator_filter_value,
    ComparisonOperators.NOT_IN: _process_membership_operator_filter_value,
    ComparisonOperators.GREATER_THAN: _process_relational_operator_filter_value,
    ComparisonOperators.GREATER_THAN_OR_EQUAL: _process_relational_operator_filter_value,
    ComparisonOperators.LESS_THAN: _process_relational_operator_filter_value,
    ComparisonOperators.LESS_THAN_OR_EQUAL: _process_relational_operator_filter_value,
}

SQL_OPERATORS_MAPPING = {
    ComparisonOperators.EQUAL: "=",
    ComparisonOperators.NOT_EQUAL: "!=",
    ComparisonOperators.IN: "in",
    ComparisonOperators.NOT_IN: "not in",
    ComparisonOperators.GREATER_THAN: ">",
    ComparisonOperators.GREATER_THAN_OR_EQUAL: ">=",
    ComparisonOperators.LESS_THAN: "<",
    ComparisonOperators.LESS_THAN_OR_EQUAL: "<=",
}


def _normalize_column_names(column_names: t.List[str]) -> t.List[str]:
    """
    Convert all column names to lowercase

    :param column_names: The column names to process.
    :return: Set of normalized columns
    """
    return [column_name.lower() for column_name in column_names]


def _remove_duplicates(column_names: t.List[str]) -> t.List[str]:
    """
    Remove duplicates from a list of column names. Case-sensitive, so it is recommended to normalize the column names
    list before using it

    :param column_names: The column names to process.
    :return: List of unique column names
    """
    seen = set()
    unique_items = []

    for s in column_names:
        if s not in seen:
            seen.add(s)
            unique_items.append(s)

    return unique_items


def _process_column_names(column_names: t.List[str]) -> t.List[str]:
    """
    Normalize column names, get rid of duplicates and validate the output list

    :param column_names: List of columns to process
    :return: Normalized and validated list of columns
    :raises ValueError: If any column does not pass the validation rule
    """
    column_names_normalized = _normalize_column_names(column_names=column_names)
    column_names_normalized_unique = _remove_duplicates(column_names_normalized)

    for column_name in column_names_normalized_unique:
        validate_column_name(column_name=column_name)

    return column_names_normalized_unique


def select(column_names: t.List[str]) -> str:
    """
    Construct a SQL SELECT body from the provided list of columns.

    :param column_names: A list of column names to be included in the SELECT clause. If the list is empty, all columns
        are selected.
    :raises ValueError: If any column name contains more than one period or is an empty string.

    :return: A string with comma-separated values of columns for subsequent use in an SQL query.

    **Example**

    Using in Mako template::
        <%!
            from casp.datastore_manager.sql import mako_helpers
        %>
        ${mako_helpers.select(["c.CELL_TYPE", "i.DATASET_ID", "c.CELL_TYPE"])}
        from `${project}.${dataset}.table_name`
        ...

    Results in::
        select c.cell_type, i.dataset_id, c.cell_type
        from `project-name-1.test_dataset.table_name`
        ...
    """
    processed_columns = _process_column_names(column_names=column_names)
    return f"select {', '.join(processed_columns)}" if processed_columns else "*"


def where(filters: t.Optional[t.Dict[str, t.Any]]) -> str:
    """
    Construct a SQL where clause from the provided filters.

    :param filters: A dictionary containing filter criteria, structured as {column_name__filter_type: value}. |br|
        Supported filter_types: |br|
            ``"eq"`` - Used for an 'equals' comparison. |br|
                Example: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``. |br|
            ``"in"`` - Used for an 'in' comparison with a set of values. |br|
                Example: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type in ('T cell', 'neuron')``.
            ``"not_eq"`` - Used for an 'not equals' comparison. Meaning that the query would exclude rows with such
                a value |br|
                Example: ``{"assay__not_eq": "Drop-seq"}
            ``"not_in"`` - Used for 'not in' comparison with a set of values to exclude. |br|
                Example: ``{"assay__not_eq": ["Drop-seq", "microwell-seq", "BD Rhapsody Targeted mRNA"]}
            ``"gt"`` - Used for an `greater than` comparison. |br|
                Example: ``{"total_mrna_umis__gt": 13000}
            ``"gte"`` - Used for an `greater than or equal` comparison. |br|
                Example: ``{"total_mrna_umis__gte": 13000}
            ``"lt"`` - Used for an `less than` comparison. |br|
                Example: ``{"total_mrna_umis__lt": 13000}
            ``"lte"`` - Used for an `less than or equal` comparison. |br|
                Example: ``{"total_mrna_umis__lte": 13000}

    :raises ValueError: If an unsupported filter type is provided.
    :return: A string representing the constructed WHERE clause, or an empty string if no valid filters are provided.

    **Example**

    Using in Mako template::
        <%!
            from casp.datastore_manager.sql import mako_helpers
        %>
        ...
        ${mako_helpers.where({"id__in": [1, 3, 5, 10], "organism__eq": "Homo sapiens"})}
        ...

    Results in::
        ...
        where
            id in (1, 3, 5, 10)
            and organism = 'Homo sapiens'
        ...
    """
    if not filters or filters is None:
        return ""

    where_conditions = []
    for filter_name, filter_value in filters.items():
        column_name, filter_type = filter_name.split("__")

        filter_value_processor = FILTER_VALUE_PROCESS_HANDLERS_MAPPING.get(filter_type)

        if not filter_value_processor:
            raise exceptions.SQLSyntaxParseException(
                f"At the moment only {', '.join(ComparisonOperators.CURRENTLY_SUPPORTED)} operators supported"
            )

        filter_value = filter_value_processor(filter_value=filter_value, filter_type=filter_type)
        sql_operator = SQL_OPERATORS_MAPPING[filter_type]

        condition = f"{column_name} {sql_operator} {filter_value}"
        where_conditions.append(condition)

    where_clause_body = "\n    and ".join(where_conditions)

    return f"where\n    {where_clause_body}" if where_clause_body else ""


def add_cell_info_required_columns(column_names: t.List[str]) -> t.List[str]:
    """
    Add required columns to the provided list of columns. Recommended to use before :func:`_process_column_names` to
    make sure all the columns are valid and unique

    :param column_names: List of columns to process
    :return: New list of column names
    """
    return [*CAS_CELL_INFO_REQUIRED_COLUMN_NAMES, *column_names]


def remove_leading_alias(column_names: t.List[str]) -> t.List[str]:
    """
    Process a list of column names by removing leading aliases.

    :param column_names: A list of column names, potentially prefixed with table aliases.
    :return: A list of column names with leading table aliases removed.

     **Example**

    >>> remove_leading_alias(["table1.column1", "alias2.column2", "column3"])
    ['column1', 'column2', 'column3']
    """
    processed_column_names = []

    for column_name in column_names:
        processed_column_names.append(column_name.split(".")[-1])

    return processed_column_names
