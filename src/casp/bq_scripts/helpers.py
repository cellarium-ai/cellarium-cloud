import typing as t

from casp.bq_scripts import constants


def _validate_column_name(column_name: str) -> None:
    """
    Validate column name. Check for number of periods in column name as well as for an empty string

    :param column_name: The column name to validate.

    :raises ValueError: If the column name contains more than one period or is an empty string.
    """
    column_split = column_name.split(".")

    if len(column_split) > 2:
        raise ValueError(
            f"Column {column_name} has more than one period in its name. It can contain only one period, "
            f"which separates the alias table name from the column name itself."
        )

    if not column_name:
        raise ValueError("Empty strings are not allowed for column names")


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


def prepare_column_names_for_extract_sql(column_names: t.Optional[t.List[str]]) -> str:
    """
    Process a list of column names by validating, normalizing, and augmenting it with required columns.
    Return a string of comma-separated, unique column names for use in SQL queries.
    Validation ensures no column name is empty and contains more than one period.

    :param column_names: The column names to process.

    :raises ValueError: If any column name contains more than one period or is an empty string.

    :return: A string with comma-separated values of columns for subsequent use in an SQL query.

    Example:
    --------
    >>> prepare_column_names_for_extract_sql(column_names=["c.CELL_TYPE", "i.DATASET_ID"])
    ["c.cas_cell_index", "c.cas_ingest_id", "c.cell_type", "i.dataset_id"]
    """
    if column_names is None or not column_names:
        return ", ".join(constants.CAS_CELL_INFO_REQUIRED_COLUMNS)

    column_names_with_required = [*constants.CAS_CELL_INFO_REQUIRED_COLUMNS, *column_names]
    column_names_normalized = _normalize_column_names(column_names=column_names_with_required)
    column_names_normalized_unique = _remove_duplicates(column_names_normalized)

    for column_name in column_names_normalized_unique:
        _validate_column_name(column_name=column_name)

    return ", ".join(column_names_normalized_unique)


def remove_leading_alias_names_extract_columns(column_names: t.List[str]) -> t.List[str]:
    """
    Process a list of column names by removing leading aliases and normalizing the names.

    This function normalizes the list of input column names, then validates each column name and removes any
    leading table alias (i.e., the "alias." prefix). It's primarily used for processing column names to be consistent
    and without aliases.

    :param column_names: A list of column names, potentially prefixed with table aliases.
    :return: A list of column names with leading table aliases removed.

    :raises ValueError: If any column name fails validation.

    Example:
    --------
    >>> remove_leading_alias_names_extract_columns(["table1.column1", "alias2.column2", "column3"])
    ['column1', 'column2', 'column3']
    """
    column_names_normalized = _normalize_column_names(column_names=column_names)
    processed_column_names = []

    for column_name in column_names_normalized:
        _validate_column_name(column_name=column_name)
        processed_column_names.append(column_name.split(".")[-1])

    return processed_column_names
