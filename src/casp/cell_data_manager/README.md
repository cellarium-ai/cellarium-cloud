# Cell Data Manager Module

This module is designed for managing cell data. Currently, we store cell data in BigQuery, which uses its own SQL dialect to handle data in the warehouse.
This module should not be limited to only BigQuery, as there is potential to switch to other data warehouse products.

For more information about BigQuery SQL dialect, please visit [this link](https://cloud.google.com/bigquery/docs/introduction-sql).

## SQL Query Management

Template generator is used for managing SQL queries due to the following reasons:

1. `cellarium-cloud` relies on a large number of queries with dynamic content (e.g., optional joins, different filters in `where` clauses).
2. To maintain the DRY (Don't Repeat Yourself) principle in SQL queries.
3. To separate SQL code from Python files, which makes our codebase cleaner, more modular, and thus more easily maintainable.

### Template Generator

[Mako template generator](https://www.makotemplates.org/) is used for this purpose. I considered both [Mako](https://www.makotemplates.org/) and [Jinja2](https://jinja.palletsprojects.com/en/2.10.x/), which were originally designed for HTML file templating. However, they seem suitable for managing SQL queries as well. You can learn more from these references:

- [Medium: Jinja + SQL = <3](https://medium.com/p/7e4dff8d8778)
- [Medium: Jinja the SQL way of Ninja](https://medium.com/analytics-and-data/jinja-the-sql-way-of-the-ninja-9a64fc815564)
- [Pushmetrics: Why Jinja and SQL?](https://pushmetrics.io/learn/jinja/why-jinja-and-sql/)

There are also wrappers around both Jinja and Mako that are specifically used for SQL management:

- [Quma (Mako wrapper)](https://github.com/ebenefuenf/quma)
- [JinjaSQL (Jinja wrapper)](https://github.com/sripathikrishnan/jinjasql)

However, it appears that these libraries are not well-maintained and/or have few stars on GitHub.

### Why Mako?

In my personal opinion, [Mako template generator](https://www.makotemplates.org/) seems to have a better syntax. It is easier to use and has a more Pythonic approach. It also allows for the seamless integration of Python helper functions within templates.

## How to use `cell_data_manager` module
1. Create a SQL template. Consider the following points:
    * Keep each query in a separate template file
    * Use `.sql.mako` extension as it explicitly points that our templates are templates and not executable queries yet
    * Utilize `cell_data_manager.sql.mako_helpers`, feel free to extend helpers if needed
2. Initialize `cell_data_manager.sql.TemplateData` instance with all the variables needed for the template
3. Call `cell_data_manager.sql.render` with a SQL template path and `cell_data_manager.sql.TemplateData` instance
4. Use a query to execute
### Example
1. SQL Template example `your_path_to_templates/templated_query.sql.mako`
```genericsql
<%!
    from casp.cell_data_manager.sql import mako_helpers
%>
select ${mako_helpers.parse_column_names(select_columns)}
from `${project}.${dataset}.cell_info` ci
% if "ii.dataset_id" in select_columns:
join `${project}.${dataset}.ingest_info` ii on (ci.ingest_id = ii.id)
$ endif
limit ${limit_by}
```
2. Create `cell_data_manager.sql.TemplateData`
```python3
from casp.cell_data_manager import sql

select_columns = ["id", "cell_type"]  # This could be used by a `mako_helpers.parse_column_names` 
additional_kwargs = {"limit_by": 100}  # Variables to be parsed in a Template
template_data = sql.TemplateData(
    project="test-project",
    dataset="my-test-dataset",
    select=select_columns,
   **additional_kwargs
)
```
3. Render SQL Query
```python3
rendered_sql_query = sql.render(
    template_path="your_path_to_templates/templated_query.sql.mako",
    template_data=template_data,
)
# Execute your query...
```
## Module Structure

The Cell Data Manager `cell_data_manager` module consists of the following structure:

```
├── data_controller # (To be developed): This should be developed during `bq_scripts` refactor
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
```
## Security Considerations

### SQL Injection Awareness

When working with SQL Template building, it's important to be aware of SQL injection vulnerabilities. SQL Template building relies on string interpolation, which can be susceptible to SQL injection attacks. Current implementation of SQL query construction doesn't leverage SQL injection checks.

**Scope**: All SQL queries executed within `bq_scripts` are designed exclusively for internal use within the `cellarium-cloud` infrastructure. `bq_scripts` is a secure module and has never been exposed to user clients or APIs. Therefore, it is safe to use template building techniques for managing our internal SQL queries.

### Best Practices

To ensure the security of your application, please follow these best practices:

1. **Avoid External Exposure**: Do not use SQL template building tools for constructing queries that will be exposed to end-users or external entities.

2. **Parameterized Queries**: If there is a need to expose SQL queries to end-users, always use [Parameterized Queries](https://cloud.google.com/bigquery/docs/parameterized-queries#python). Parameterized queries help prevent SQL injection by separating data from SQL code.

3. **Consider extending `mako_helpers`**: To simplify the implementation of parameterized queries, consider expanding `mako_helpers` to work seamlessly with [Parameterized Queries](https://cloud.google.com/bigquery/docs/parameterized-queries#python).

By following these practices, you can ensure the safe and secure use of this module within the project.

## To Consider

Currently, there's no SQL formatting in the `casp.cell_data_manager.sql` module. There's no SQL formatter tool that supports the Mako template manager. However, there is one that supports Jinja:
(sqlfmt)[https://github.com/tconbeer/sqlfmt]. I didn't want to use Jinja just because there's a formatting tool that supports it. Hopefully, we will find a SQL formatting tool suitable for our use case, or we might consider contributing to this repository to write a Mako formatting connector.
At the moment, it is advised to keep all `sql.mako` templates in a similar format manually.