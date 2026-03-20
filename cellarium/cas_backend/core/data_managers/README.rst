Data Managers
=============

The ``cellarium.cas_backend.core.data_managers`` package contains warehouse-facing data access code used by compute
workflows. It encapsulates SQL template rendering, BigQuery-oriented query helpers, and module-specific data access
objects.

What Lives Here
---------------

- ``base_data_manager.py``: shared base class for data managers
- ``cell_operations.py``: data access for compute workflows
- ``cell_quota.py``: quota-related access helpers
- ``cellarium_general.py``: shared metadata and schema lookups
- ``sql/``: SQL rendering, template data, helpers, and validation logic
- ``sql_templates/`` and ``sql_queries/``: query templates used by the data managers

Public Entry Points
-------------------

- ``cellarium.cas_backend.core.data_managers.base_data_manager.BaseDataManager``
- ``cellarium.cas_backend.core.data_managers.sql.render``
- ``cellarium.cas_backend.core.data_managers.sql.TemplateData``

Dependencies
------------

Data managers are used primarily by the compute service and depend on:

- BigQuery and related Google Cloud clients
- shared configuration from ``cellarium.cas_backend.core.config``
- optional database session access from the core db package

SQL Template Approach
---------------------

Dynamic SQL is stored in ``.sql.mako`` templates so query logic stays close to the data layer and out of the app
handlers. This is useful for optional joins, conditional filters, and composable warehouse queries.

When adding a new query:

1. Create a dedicated ``.sql.mako`` template file.
2. Pass all required values through ``TemplateData``.
3. Render with ``sql.render(...)`` from the data manager layer.
4. Execute the rendered query through the appropriate client.

Security Notes
--------------

These templates are intended for internal query generation. Do not use raw string interpolation patterns here for
externally supplied ad hoc SQL. If a workflow needs end-user parameters, keep them constrained and use the BigQuery
parameterization features where applicable.

Related Documentation
---------------------

- ``cellarium/cas_backend/core/db/README.rst``
- ``docs/source/modules/running_locally.rst``
