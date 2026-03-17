class BaseSQLTemplateException(Exception):
    pass


class SQLSyntaxParseException(BaseSQLTemplateException):
    pass


class UnsupportedSQLTypeException(BaseSQLTemplateException):
    pass
