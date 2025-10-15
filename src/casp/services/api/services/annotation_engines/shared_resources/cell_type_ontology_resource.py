import json
import typing as t

from smart_open import open

from casp.services import settings


class CellOntologyResource:
    """
    Handles cell ontology resource data, such as the ancestors dictionary and cell ontology term ID mappings.

    Attributes:
        ancestors_dictionary: Maps cell ontology term IDs to their ancestor IDs.
        ontology_term_id_to_name_dict: Maps cell ontology term IDs to cell type names.

    Raises:
        ValueError: If either ``ancestors_dictionary`` or ``cell_ontology_term_id_to_cell_type`` is missing from the
            provided dictionary.
    """

    def __init__(self, cell_ontology_resource_dict: t.Optional[t.Dict[str, t.Any]] = None):
        if cell_ontology_resource_dict is None:
            with open(settings.GCS_CELL_ONTOLOGY_RESOURCE_FILE, "rb") as f:
                cell_ontology_resource_dict = json.loads(f.read())

        if "ancestors_dictionary" not in cell_ontology_resource_dict:
            raise ValueError("`ancestors_dictionary` is not found in the cell ontology resource dictionary")
        if "cell_ontology_term_id_to_cell_type" not in cell_ontology_resource_dict:
            raise ValueError(
                "`cell_ontology_term_id_to_cell_type` mapping is not found in the cell ontology resource dictionary"
            )

        self.ancestors_dictionary = cell_ontology_resource_dict["ancestors_dictionary"]
        self.ontology_term_id_to_name_dict = cell_ontology_resource_dict["cell_ontology_term_id_to_cell_type"]
