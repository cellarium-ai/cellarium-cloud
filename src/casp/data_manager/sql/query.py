import logging
import typing as t

import sqlparse
from mako.template import Template

from casp.data_manager.sql.template_data import TemplateData
from casp.data_manager.sql.validation.query_validators import BigQuerySQLSyntaxValidator

if t.TYPE_CHECKING:
    from casp.data_manager.sql.validation.query_validators.base_query_validator import SQLSyntaxValidator

logger = logging.getLogger(__name__)


def render(
    template_path: str,
    template_data: TemplateData,
    sql_query_validator_class: "SQLSyntaxValidator" = BigQuerySQLSyntaxValidator,
    sql_query_validator_on: bool = False,
) -> str:
    """
    Render a SQL query from a template file, populate it with provided data, and optionally validate the generated
    query's syntax.

    This function is designed to read a specified SQL template file, populate it with the data encapsulated in a
    :class:`~datastore_manager.sql.TemplateData` instance, and perform an optional syntax validation on the resulting
    SQL query using the designated validator class. Enabling syntax validation is recommended during development to
    identify potential syntax issues proactively, but not recommended in prod mode as it might cost additional
    computations and additional database hits.

    :param template_path: The file path to the SQL template.
    :param template_data: A :class:`~datastore_manager.sql.TemplateData` instance with data to be used in the SQL
        template rendering.
    :param sql_query_validator_class: SQL query syntax validator class. Must be a subclass
        of :class:`~datastore_manager.sql.validation.query_validators.SQLSyntaxValidator` |br|
        `Default:` :class:`~datastore_manager.sql.validation.query_validators.BigQuerySQLSyntaxValidator`
    :param sql_query_validator_on: A flag describing whether to use a sql query validation. Recommended to use only
        in dev mode |br|
        `Default:` ``False``
    :return: A string representing the rendered SQL query.
    """
    template = Template(filename=template_path)
    rendered_sql_query = template.render(**template_data.data)
    rendered_sql_query = rendered_sql_query.lstrip()
    rendered_sql_query = sqlparse.format(rendered_sql_query, reindent=True, keyword_case="lower")
    if sql_query_validator_on:
        sql_query_validator_class.validate_syntax(sql_query=rendered_sql_query)

    logger.info(f"Rendered SQL Query:\n{rendered_sql_query}")
    return rendered_sql_query
