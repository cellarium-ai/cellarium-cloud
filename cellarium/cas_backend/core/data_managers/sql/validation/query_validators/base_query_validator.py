class SQLSyntaxValidator:
    """
    Base abstract class for SQL query validation. This class should not be instantiated directly.

    Subclasses are required to implement the `validate_syntax` method to provide specific validation logic.
    """

    @classmethod
    def validate_syntax(cls, sql_query: str) -> None:
        """
        Validate query syntax and raises an error if the query is not valid.

        Implementers should: Raise a `exceptions.ValidationError` if the query does not pass validation.

        :param sql_query: The SQL query to validate as a string.
        :raises exceptions.ValidationError: If the query does not pass validation.
        """
        raise NotImplementedError
