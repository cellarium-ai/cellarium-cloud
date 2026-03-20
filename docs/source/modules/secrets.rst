Secret Variables
================

.. _Project Secrets:

Project Runtime Secrets
-----------------------

Runtime configuration is loaded from ``settings/.env`` in local development and mounted into containers at
``/app/settings/.env`` in Cloud Run.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Secret Name
     - Description
   * - ``ENVIRONMENT``
     - Runtime environment such as ``local``, ``test``, or ``production``
   * - ``GOOGLE_SERVICE_ACCOUNT_CREDENTIALS``
     - Optional JSON credentials blob for Google Cloud clients
   * - ``PROJECT_BUCKET_NAME``
     - Bucket used for shared project resources
   * - ``DB_HOST``
     - Local PostgreSQL host
   * - ``DB_PORT``
     - PostgreSQL port
   * - ``DB_NAME``
     - Database name
   * - ``DB_USER``
     - Database username
   * - ``DB_PASSWORD``
     - Database password
   * - ``DB_PRIVATE_IP``
     - Private IP for Cloud SQL connectivity
   * - ``DB_CONNECTION_NAME``
     - Cloud SQL connection name
   * - ``DB_INSTANCE_UNIX_SOCKET``
     - Unix socket path for Cloud SQL
   * - ``FLASK_SECRET_KEY``
     - Flask session secret for admin
   * - ``FLASK_SECURITY_PASSWORD_SALT``
     - Password salt for admin
   * - ``FLASK_BASIC_AUTH_USERNAME``
     - Basic-auth username for admin
   * - ``FLASK_BASIC_AUTH_PASSWORD``
     - Basic-auth password for admin
   * - ``SENDGRID_API_KEY``
     - SendGrid API key
   * - ``FROM_ADDRESS``
     - Default sender address
   * - ``SENTRY_DSN``
     - Sentry DSN
   * - ``SENTRY_ENVIRONMENT``
     - Sentry environment tag
   * - ``SENTRY_ENABLE_TRACING``
     - Enable request tracing
   * - ``SENTRY_TRACES_SAMPLE_RATE``
     - Trace sample rate
   * - ``SENTRY_PROFILES_SAMPLE_RATE``
     - Profile sample rate

GitHub Actions Secrets
----------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Secret Name
     - Description
   * - ``GCP_PROJECT_ID``
     - Google Cloud project ID
   * - ``GCP_PROVIDER_ID``
     - Workload identity provider ID
   * - ``GCP_SERVICE_ACCOUNT_EMAIL``
     - Deployment service account email
   * - ``CAS_SERVICES_ENV_B64_DEV``
     - Base64-encoded environment file for development deployments
   * - ``CAS_SERVICES_ENV_B64_PROD``
     - Base64-encoded environment file for production deployments
   * - ``TEST_API_KEY``
     - API key used by integration-test automation

Notes
-----

The authoritative list of runtime fields is ``cellarium.cas_backend.core.config``. Update this page when new required
environment variables are introduced for developers or CI/CD.
