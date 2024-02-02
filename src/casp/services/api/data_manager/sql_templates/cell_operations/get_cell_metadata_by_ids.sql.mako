<%
"""
Description:
  This template is used for generating a SQL query for accessing the 'cas_cell_info' table.
  It utilizes helper functions from 'casp.data_manager.sql.mako_helpers' for selecting
  columns and applying filters.

Modules and Functions Used:
  - mako_helpers (mh): A module providing helper functions like 'select' and 'where' for
    constructing SQL queries.

Parameters:
  - select_columns: A list or tuple of column names to be selected in the SQL query.
  - filter_statements: Conditions for filtering the data. Please refer to `mako_helpers.where`
    for more details on the filter statement format.

Tables Accessed:
  - cas_cell_info: A table within the specified project and dataset in BigQuery.

Output Format:
  It would depend on what columns are specified in the 'select_columns' parameter.

Notes:
  - The 'select' function from mako_helpers is used to specify the columns to be retrieved.
  - The 'where' function applies the given filter conditions to the query.
"""
%>

<%
    from casp.data_manager.sql import mako_helpers as mh
%>
${mh.select(select_columns)}
from `${project}.${dataset}.cas_cell_info` ci
${mh.where(filter_statements)}
