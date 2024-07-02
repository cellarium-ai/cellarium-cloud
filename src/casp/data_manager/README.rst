Data Manager
============

..
    Are we still using bigquery?  Will we still be using it at the time when we start letting people actually use this thing?

This module is designed for managing single-cell data in Cellarium DataStore. Currently, single-cell data are stored in BigQuery, which uses its own SQL dialect to handle data in the warehouse. This module should not be limited to only BigQuery, as there is potential to switch to other data warehouse products.

For more information about BigQuery SQL dialect, please visit `BigQuery SQL dialect <https://cloud.google.com/bigquery/docs/introduction-sql>`_.

SQL Query Management
--------------------

The template generator is used for managing SQL queries due to the following reasons:

1. ``cellarium-cloud`` relies on a large number of queries with dynamic content (e.g., optional joins, different filters in ``where`` clauses).
2. To maintain the DRY (Don't Repeat Yourself) principle in SQL queries.
3. To separate SQL code from Python files, which makes our codebase cleaner, more modular, and thus more easily maintainable.

Template Generator
~~~~~~~~~~~~~~~~~~

`Mako template generator <https://www.makotemplates.org/>`_ is chosen for this purpose. Both `Mako <https://www.makotemplates.org/>`_ and `Jinja2 <https://jinja.palletsprojects.com/en/2.10.x/>`_, originally designed for HTML file templating, were considered suitable for managing SQL queries. More information can be found in these references:

- `Medium: Jinja + SQL = <3 <https://medium.com/p/7e4dff8d8778>`_
- `Medium: Jinja the SQL way of Ninja <https://medium.com/analytics-and-data/jinja-the-sql-way-of-the-ninja-9a64fc815564>`_
- `Pushmetrics: Why Jinja and SQL? <https://pushmetrics.io/learn/jinja/why-jinja-and-sql/>`_

Wrapper libraries for both Jinja and Mako specifically used for SQL management are also available:

- `Quma (Mako wrapper) <https://github.com/ebenefuenf/quma>`_
- `JinjaSQL (Jinja wrapper) <https://github.com/sripathikrishnan/jinjasql>`_

However, it appears that these libraries are not well-maintained and/or have few stars on GitHub.

Why Mako?
~~~~~~~~~

The `Mako template generator <https://www.makotemplates.org/>`_ is preferred for its better syntax. It is found to be easier to use and has a more Pythonic approach. It also allows for the seamless integration of Python helper functions within templates.

Using the `data_access_manager` module
--------------------------------------

1. Create a SQL template. Consider the following points:
   * Keep each query in a separate template file.
   * Use ``.sql.mako`` extension as it explicitly indicates that the templates are not executable queries yet.
   * Utilize ``data_access_manager.sql.mako_helpers``, feel free to extend helpers if needed.

2. Initialize ``data_access_manager.sql.TemplateData`` instance with all the variables needed for the template.

3. Call ``data_access_manager.sql.render`` with a SQL template path and ``data_access_manager.sql.TemplateData`` instance.

4. Use a query to execute.

Example
~~~~~~~

1. SQL Template example ``your_path_to_templates/templated_query.sql.mako``

   .. code-block::

       <%!
           from casp.data_access_manager.sql import mako_helpers
       %>
       select ${mako_helpers.parse_column_names(select_columns)}
       from `${project}.${dataset}.cell_info` ci
       % if "ii.dataset_id" in select_columns:
       join `${project}.${dataset}.ingest_info` ii on (ci.ingest_id = ii.id)
       $ endif
       limit ${limit_by}

2. Create ``data_access_manager.sql.TemplateData``

   .. code-block:: python

       from casp.data_manager import sql

       select_columns = ["id", "cell_type"]  # This could be used by a `mako_helpers.parse_column_names` 
       additional_kwargs = {"limit_by": 100}  # Variables to be parsed in a Template
       template_data = sql.TemplateData(
           project="test-project",
           dataset="my-test-dataset",
           select=select_columns,
           **additional_kwargs
       )

3. Render SQL Query

   .. code-block:: python

       rendered_sql_query = sql.render(
           template_path="your_path_to_templates/templated_query.sql.mako",
           template_data=template_data,
       )
       # Execute your query...

Module Structure
----------------

The Data Manager ``data_manager`` module consists of the following structure:

.. code-block:: text

    ├── base_data_manager.py # Base class for all data managers in the project
    └── sql
        ├── __init__.py
        ├── constants.py
        ├── mako_helpers.py  # Helper functions used in SQL templates
        ├── query.py  # Query management / validation directory
        ├── template_data.py # Template Data Parsing class
        └── templates # All SQL templates go here
            └── ... 
        └── validation
            └── ...                

Security Considerations
-----------------------

SQL Injection Awareness
~~~~~~~~~~~~~~~~~~~~~~~

When working with SQL Template building, it's important to be aware of SQL injection vulnerabilities. SQL Template building relies on string interpolation, which can be susceptible to SQL injection attacks. Current implementation of SQL query construction doesn't leverage SQL injection checks.

**Scope**: All SQL queries executed within ``bq_scripts`` are designed exclusively for internal use within the ``cellarium-cloud`` infrastructure. ``bq_scripts`` is a secure module and has never been exposed to user clients or APIs. Therefore, it is safe to use template building techniques for managing our internal SQL queries.

Best Practices
~~~~~~~~~~~~~~

To ensure the security of your application, please follow these best practices:

1. **Avoid External Exposure**: Do not use SQL template building tools for constructing queries that will be exposed to end-users or external entities.

2. **Parameterized Queries**: If there is a need to expose SQL queries to end-users, always use `Parameterized Queries <https://cloud.google.com/bigquery/docs/parameterized-queries#python>`_. Parameterized queries help prevent SQL injection by separating data from SQL code.

3. **Consider extending `mako_helpers`**: To simplify the implementation of parameterized queries, consider expanding `mako_helpers` to work seamlessly with `Parameterized Queries <https://cloud.google.com/bigquery/docs/parameterized-queries#python>`_.

By following these practices, you can ensure the safe and secure use of this module within the project.

To Consider
-----------

Currently, there is no SQL formatting in the ``casp.data_access_manager.sql`` module. There's no SQL formatter tool that supports the Mako template manager. However, there is one that supports Jinja: `sqlfmt <https://github.com/tconbeer/sqlfmt>`_. Jinja was not chosen solely for the reason that a formatting tool supports it. Hopefully, a SQL formatting tool suitable for SQL Mako templates will be found. At the moment, it is advised to keep all ``sql.mako`` templates in a similar format manually. Please use lower-case style for SQL query consistency.
