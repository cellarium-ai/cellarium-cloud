class ValidationError(Exception):
    pass


class QueryValidationError(ValidationError):
    pass


class TemplateDataValidationError(ValidationError):
    pass
