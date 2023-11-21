import typing as t

from google.cloud import bigquery
from google.cloud import exceptions as gc_exceptions

from casp.datastore_manager.sql.validation import exceptions
from casp.datastore_manager.sql.validation.query_validators.base_query_validator import SQLSyntaxValidator

if t.TYPE_CHECKING:
    import google


class BigQuerySQLSyntaxValidator(SQLSyntaxValidator):
    """
    Validate BigQuery SQL syntax using the dry-run feature provided by BigQuery.

    The dry-run feature allows the validation of SQL syntax by the BigQuery backend without actually executing
    the query, thus avoiding query execution costs and side effects.

    For more details on the dry-run feature, refer to:
    https://cloud.google.com/bigquery/docs/running-queries#dry-run
    """

    @staticmethod
    def _bq_dry_run_query(query) -> "google.cloud.bigquery.job.QueryJob":
        """
        Perform a BigQuery Dry Run query to validate the syntax of a SQL query.

        :param query: The SQL query to be validated.
        :returns: A BigQuery QueryJob object representing the validation job.

        :raises google.api_core.exceptions.GoogleCloudError: If there is an error during the Dry Run query execution,
            typically due to syntax errors or BigQuery service issues.
        """
        job_config = bigquery.QueryJobConfig(dry_run=True, use_query_cache=True)

        bq_client = bigquery.Client()
        return bq_client.query(
            query,
            job_config=job_config,
        )

    @staticmethod
    def _parse_google_bq_error_message(error: "google.cloud.exceptions.GoogleCloudError") -> str:
        """
        Extract the meaningful error message from a GoogleCloudError exception.

        This method parses the error message to extract the core information relevant to the SQL syntax error
        encountered during the dry-run.

        :param error: The :class:`google.cloud.exceptions.GoogleCloudError` exception captured during the dry-run
            query execution.
        :return: A string representing the simplified error message relevant for SQL syntax validation.
        """
        error_minus_job_url = str(error).split("jobs?prettyPrint=false:")
        split_error = error_minus_job_url[1].split("\n\n")
        return split_error[0]

    @classmethod
    def validate_syntax(cls, sql_query: str) -> None:
        """
        Validate the syntax of a rendered SQL query using BigQuery Dry Run. Catch
        :class:`google.cloud.exceptions.GoogleCloudError` error and raise :class:`exceptions.QueryValidationError` with
        a parsed and simplified syntax error message from Google Cloud API.

        :param sql_query: Rendered SQL query to be validated.

        :raises exceptions.ValidationError: If there is a syntax error in the query.
        """
        try:
            cls._bq_dry_run_query(sql_query)
        except gc_exceptions.GoogleCloudError as e:
            syntax_error_message = cls._parse_google_bq_error_message(e)
            raise exceptions.QueryValidationError(syntax_error_message)
