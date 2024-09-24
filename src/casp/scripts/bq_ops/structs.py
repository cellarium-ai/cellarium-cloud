import typing as t
from dataclasses import dataclass


@dataclass
class CategoricalColumnValues:
    column_name: str
    unique_values: t.List[str]
