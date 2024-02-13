import typing as t

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


def _int_value_processor(v: int) -> str:
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
    Construct a SQL WHERE clause from the provided filters.

    :param filters: A dictionary containing filter criteria, structured as {column_name__filter_type: value}. |br|
        Supported filter_types: |br|
            ``"eq"`` - Used for an 'equals' comparison. |br|
                Example: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``. |br|
            ``"in"`` - Used for an 'in' comparison with a set of values. |br|
                Example: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type in ('T cell', 'neuron')``.
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

        if filter_type == ComparisonOperators.EQUAL:
            if isinstance(filter_value, str):
                filter_value = _string_value_processor(v=filter_value)
            if isinstance(filter_value, bool):
                filter_value = _bool_value_processor(v=filter_value)

            condition = f"{column_name} = {filter_value}"

        elif filter_type == ComparisonOperators.IN:
            if isinstance(filter_value[0], str):
                filter_value = [_string_value_processor(x) for x in filter_value]
            if isinstance(filter_value[0], int):
                filter_value = [_int_value_processor(x) for x in filter_value]

            filter_value = ", ".join(filter_value)
            condition = f"{column_name} in ({filter_value})"
        else:
            raise ValueError(
                f"At the moment only {', '.join(ComparisonOperators.CURRENTLY_SUPPORTED)} operators supported"
            )

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
