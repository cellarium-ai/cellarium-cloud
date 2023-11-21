import typing as t

from casp.datastore_manager.sql.constants import TemplateDataDictNames
from casp.datastore_manager.sql.validation import template_data_validator


class TemplateData:
    """
    Represents data for SQL templates.

    This class provides a structured way to store and pass data for SQL templates.

    :param project: The BigQuery project.
    :param dataset: The BigQuery dataset.
    :param select: The SELECT statement. |br|
        `Default:` ``None``
    :param filters: A dictionary containing filter criteria, structured as {column_name__filter_type: value}. |br|
        Supported filter_types: |br|
            ``"eq"`` - Used for an 'equals' comparison. |br|
                Example: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``. |br|
            ``"in"`` - Used for an 'in' comparison with a set of values. |br|
                Example: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type in
                    ('T cell', 'neuron')``. |br|
        `Default:` ``None``
    :param other_kwargs: Additional keyword arguments to include in the data.

    :raises ValueError: If project or dataset name is not specified.

    Example
    -------
    Create a TemplateData instance with project, dataset and filters:

    >>> template_data = TemplateData(
    >>>     project="your-project",
    >>>     dataset="your_dataset",
    >>>     filters={"c.organism__eq": "Homo sapiens", "c.cell_type__in": ["T Cell", "B Cell", "neuron"]}
    >>> )
    """

    def __init__(
        self,
        project: str,
        dataset: str,
        select: t.Optional[t.List[str]] = None,
        filters: t.Optional[t.Dict[str, t.Any]] = None,
        **other_kwargs: t.Union[str, bool, int, float],
    ) -> None:
        template_data_validator.validate_not_empty(value=project)
        template_data_validator.validate_not_empty(value=dataset)

        self.data = {
            TemplateDataDictNames.PROJECT: project,
            TemplateDataDictNames.DATASET: dataset,
        }

        if select is not None:
            template_data_validator.validate_not_empty(value=select)
            self.data[TemplateDataDictNames.SELECT] = select
        if filters is not None:
            template_data_validator.validate_filters(filters=filters)
            self.data[TemplateDataDictNames.FILTERS] = filters

        for k, v in other_kwargs.items():
            self.data[k] = v
