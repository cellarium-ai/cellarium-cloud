Secret Variables
=================

.. _Project Secrets:

Project Secrets
---------------
The following secrets are required for the project to build and deploy. They have to be located in the
`settings/.env` file. The file is not included in the repository for security reasons. The file should be
created manually.

+------------------------------------+-----------------------------------------------------------------------+
| Secret Name                        | Description                                                           |
+====================================+=======================================================================+
| GOOGLE_SERVICE_ACCOUNT_CREDENTIALS | (deprecated) A string dump of credentials for the service account used|
|                                    | to access Google Cloud Storage.                                       |
+------------------------------------+-----------------------------------------------------------------------+
| NEPTUNE_API_KEY                    | API key for Neptune (Needed only for ML pipelines)                    |
+------------------------------------+-----------------------------------------------------------------------+
| DB_PORT                            | Port for the database                                                 |
+------------------------------------+-----------------------------------------------------------------------+
| DB_NAME                            | Name of the database                                                  |
+------------------------------------+-----------------------------------------------------------------------+
| DB_USER                            | Username for the database                                             |
+------------------------------------+-----------------------------------------------------------------------+
| DB_PASSWORD                        | Password for the database                                             |
+------------------------------------+-----------------------------------------------------------------------+
| DB_PRIVATE_IP                      | Private IP for the Cloud SQL DB                                       |
+------------------------------------+-----------------------------------------------------------------------+
| DB_CONNECTION_NAME                 | Connection name for the database (Not needed for local env)           |
+------------------------------------+-----------------------------------------------------------------------+
| DB_INSTANCE_UNIX_SOCKET            | Unix socket for the database (Not needed for local env)               |
+------------------------------------+-----------------------------------------------------------------------+
| FLASK_SECRET_KEY                   | Secret key for Flask server                                           |
+------------------------------------+-----------------------------------------------------------------------+
| FLASK_SECURITY_PASSWORD_SALT       | Password salt for Flask server                                        |
+------------------------------------+-----------------------------------------------------------------------+
| FLASK_BASIC_AUTH_USERNAME          | Username for basic auth in admin                                      |
+------------------------------------+-----------------------------------------------------------------------+
| FLASK_BASIC_AUTH_PASSWORD          | Password for basic auth in admin                                      |
+------------------------------------+-----------------------------------------------------------------------+
| ENVIRONMENT                        | Environment for the project `development`, `production` or `local`    |
+------------------------------------+-----------------------------------------------------------------------+
| PROJECT_BUCKET_NAME                | Name of the bucket for the project                                    |
+------------------------------------+-----------------------------------------------------------------------+
| SENDGRID_API_KEY                   | Key to use to send emails using the Sendgrid API                      |
+------------------------------------+-----------------------------------------------------------------------+
| FROM_ADDRESS                       | Email address to send emails from                                     |
+------------------------------------+-----------------------------------------------------------------------+
| SENTRY_DSN                         | Connection string for reporting to Sentry. Must be set to use Sentry  |
+------------------------------------+-----------------------------------------------------------------------+
| SENTRY_ENVIRONMENT                 | Environment that Sentry events will be associated with.  Must be set |
|                                    | to use Sentry                                                         |
+------------------------------------+-----------------------------------------------------------------------+
| SENTRY_ENABLE_TRACING              | If true, transactions and trace data will be generated and captured   |
+------------------------------------+-----------------------------------------------------------------------+
| SENTRY_TRACES_SAMPLE_RATE          | The percentage chance a given transaction will be sent to Sentry      |
+------------------------------------+-----------------------------------------------------------------------+
| SENTRY_PROFILES_SAMPLE_RATE        | Sample rate in Sentry relative to SENTRY_TRACES_SAMPLE_RATE           |
+------------------------------------+-----------------------------------------------------------------------+



Github Actions Required Secrets
-------------------------------
These secrets are used by GitHub Actions to connect with Google Cloud Platform and update docker images

+---------------------------+-----------------------------------------------------------------------+
| Secret Name               | Description                                                           |
+===========================+=======================================================================+
| GCP_PROJECT_ID            | Google Cloud Platform project ID                                      |
+---------------------------+-----------------------------------------------------------------------+
| GCP_PROVIDER_ID           | Google Cloud Platform provider ID                                     |
+---------------------------+-----------------------------------------------------------------------+
| GCP_SERVICE_ACCOUNT_EMAIL | Google Cloud Platform service account email                           |
+---------------------------+-----------------------------------------------------------------------+
| GCP_SERVICE_ACCOUNT_KEY   | Google Cloud Platform service account key                             |
+---------------------------+-----------------------------------------------------------------------+
| CAS_SERVICES_ENV_B64      | Base64 encoded `.env` file which is passed to docker's secret         |
|                           | environment. Should include secrets from :ref:`Project Secrets`       |
+---------------------------+-----------------------------------------------------------------------+
