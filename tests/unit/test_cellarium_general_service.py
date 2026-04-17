"""
Tests things in the cellarium_general_service that can reasonably be tested with a unit test.
"""

import json
import typing as t

from mockito import mock, unstub
import pytest

from cellarium.cas_backend.apps.compute.schemas import CellOntologyResourceResponse
from cellarium.cas_backend.apps.compute.schemas.cellarium_general import CASModel as CASModelSchema
from cellarium.cas_backend.apps.compute.schemas.cellarium_general import OntologicalColumnInfo
from cellarium.cas_backend.apps.compute.services.cellarium_general_service import CellariumGeneralService
from cellarium.cas_backend.apps.compute.services.consensus_engine.strategies.ontology_aware import CellOntologyResource
from cellarium.cas_backend.apps.compute.services.exceptions import InvalidClientVersionException
from cellarium.cas_backend.core import constants


@pytest.fixture
def cellarium_general_service():
    service = CellariumGeneralService(cellarium_general_dm=mock())
    yield service
    unstub()


@pytest.fixture
def ontology_resource_dict() -> dict[str, t.Any]:
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


@pytest.fixture
def ontology_resource(ontology_resource_dict) -> CellOntologyResource:
    return CellOntologyResource(cell_ontology_resource_dict=ontology_resource_dict)


def test_validate_client_version_valid(cellarium_general_service):
    """A version matching the minimum required version is accepted."""
    assert cellarium_general_service.validate_client_version(constants.MIN_CLIENT_VERSION)


def test_validate_client_version_too_old(cellarium_general_service):
    """A version below the minimum returns False without raising."""
    assert not cellarium_general_service.validate_client_version("1.3.0")


def test_validate_client_version_invalid(cellarium_general_service):
    """A non-semver string raises InvalidClientVersionException."""
    with pytest.raises(InvalidClientVersionException):
        cellarium_general_service.validate_client_version("invalid_version")


def test_cl_names_is_sorted(ontology_resource):
    """cl_names are in lexicographic order for stable downstream enumeration."""
    assert ontology_resource.cl_names == sorted(ontology_resource.cl_names)


def test_cl_names_contains_all_term_ids(ontology_resource, ontology_resource_dict):
    """cl_names covers exactly the same terms as cell_ontology_term_id_to_cell_type."""
    assert set(ontology_resource.cl_names) == set(ontology_resource_dict["cell_ontology_term_id_to_cell_type"].keys())


def test_children_dictionary_loaded(ontology_resource, ontology_resource_dict):
    """children_dictionary is populated and matches the raw fixture data."""
    assert ontology_resource.children_dictionary is not None
    assert ontology_resource.children_dictionary == ontology_resource_dict["children_dictionary"]


def test_shortest_path_lengths_loaded(ontology_resource):
    """shortest_path_lengths_from_cell_root is populated and root has length 0."""
    assert ontology_resource.shortest_path_lengths_from_cell_root is not None
    assert ontology_resource.shortest_path_lengths_from_cell_root["CL:0000000"] == 0


def test_longest_path_lengths_loaded(ontology_resource):
    """longest_path_lengths_from_cell_root is populated and root has length 0."""
    assert ontology_resource.longest_path_lengths_from_cell_root is not None
    assert ontology_resource.longest_path_lengths_from_cell_root["CL:0000000"] == 0


def test_new_fields_optional_when_missing(ontology_resource_dict):
    """CellOntologyResource can be built from a legacy dict missing the three new fields."""
    minimal_dict = {
        "ancestors_dictionary": ontology_resource_dict["ancestors_dictionary"],
        "cell_ontology_term_id_to_cell_type": ontology_resource_dict["cell_ontology_term_id_to_cell_type"],
    }
    resource = CellOntologyResource(cell_ontology_resource_dict=minimal_dict)
    assert resource.children_dictionary is None
    assert resource.shortest_path_lengths_from_cell_root is None
    assert resource.longest_path_lengths_from_cell_root is None


def test_response_has_correct_keys(ontology_resource):
    """Serialized response contains exactly the five public fields; ancestors_dictionary is excluded."""
    response = CellOntologyResourceResponse(
        cl_names=ontology_resource.cl_names,
        cell_ontology_term_id_to_cell_type=ontology_resource.ontology_term_id_to_name_dict,
        children_dictionary=ontology_resource.children_dictionary,
        shortest_path_lengths_from_cell_root=ontology_resource.shortest_path_lengths_from_cell_root,
        longest_path_lengths_from_cell_root=ontology_resource.longest_path_lengths_from_cell_root,
    )
    dumped = response.model_dump()
    assert set(dumped.keys()) == {
        "cl_names",
        "cell_ontology_term_id_to_cell_type",
        "children_dictionary",
        "shortest_path_lengths_from_cell_root",
        "longest_path_lengths_from_cell_root",
    }
    assert "ancestors_dictionary" not in dumped


def test_cl_names_sorted_in_response(ontology_resource):
    """cl_names in the serialized response are in lexicographic order."""
    response = CellOntologyResourceResponse(
        cl_names=ontology_resource.cl_names,
        cell_ontology_term_id_to_cell_type=ontology_resource.ontology_term_id_to_name_dict,
        children_dictionary=ontology_resource.children_dictionary,
        shortest_path_lengths_from_cell_root=ontology_resource.shortest_path_lengths_from_cell_root,
        longest_path_lengths_from_cell_root=ontology_resource.longest_path_lengths_from_cell_root,
    )
    assert response.cl_names == sorted(response.cl_names)


def test_ontological_column_info_from_attributes():
    """OntologicalColumnInfo is populated correctly from object attributes."""

    class FakeOntologicalColumn:
        ontology_resource_name = "cell_type"
        column_name = "cell_type"
        description = "Cell type ontology"

    info = OntologicalColumnInfo.model_validate(FakeOntologicalColumn(), from_attributes=True)
    assert info.ontology_resource_name == "cell_type"
    assert info.column_name == "cell_type"
    assert info.description == "Cell type ontology"


def test_ontological_column_info_description_optional():
    """OntologicalColumnInfo description defaults to None when not set."""

    class FakeOntologicalColumn:
        ontology_resource_name = "disease"
        column_name = "disease"
        description = None

    info = OntologicalColumnInfo.model_validate(FakeOntologicalColumn(), from_attributes=True)
    assert info.description is None


def test_cas_model_ontological_columns_populated():
    """CASModel reads ontological_columns via AliasPath through cell_info_metadata."""

    class FakeColumn:
        ontology_resource_name = "cell_type"
        column_name = "cell_type"
        description = "Cell type ontology"

    class FakeCellInfoMetadata:
        ontological_columns = [FakeColumn()]

    class FakeCASModel:
        model_name = "human-pca-001"
        description = "Test model"
        schema_name = "refdata-gex-GRCh38-2020-A"
        is_default_model = False
        embedding_dimension = 512
        cell_info_metadata = FakeCellInfoMetadata()

    schema = CASModelSchema.model_validate(FakeCASModel(), from_attributes=True)
    assert len(schema.ontological_columns) == 1
    assert schema.ontological_columns[0].ontology_resource_name == "cell_type"
    assert schema.ontological_columns[0].column_name == "cell_type"
    assert schema.ontological_columns[0].description == "Cell type ontology"


def test_cas_model_ontological_columns_empty_list():
    """CASModel ontological_columns defaults to an empty list when cell_info_metadata has none."""

    class FakeCellInfoMetadata:
        ontological_columns = []

    class FakeCASModel:
        model_name = "human-pca-001"
        description = "Test model"
        schema_name = "refdata-gex-GRCh38-2020-A"
        is_default_model = False
        embedding_dimension = 512
        cell_info_metadata = FakeCellInfoMetadata()

    schema = CASModelSchema.model_validate(FakeCASModel(), from_attributes=True)
    assert schema.ontological_columns == []


def test_cas_model_multiple_ontological_columns():
    """CASModel includes all ontological columns when multiple are associated."""

    class FakeColumn:
        def __init__(self, name, col, desc):
            self.ontology_resource_name = name
            self.column_name = col
            self.description = desc

    class FakeCellInfoMetadata:
        ontological_columns = [
            FakeColumn("cell_type", "cell_type", "Cell type"),
            FakeColumn("disease", "disease", None),
        ]

    class FakeCASModel:
        model_name = "human-pca-001"
        description = "Test model"
        schema_name = "refdata-gex-GRCh38-2020-A"
        is_default_model = False
        embedding_dimension = 512
        cell_info_metadata = FakeCellInfoMetadata()

    schema = CASModelSchema.model_validate(FakeCASModel(), from_attributes=True)
    assert len(schema.ontological_columns) == 2
    names = {c.ontology_resource_name for c in schema.ontological_columns}
    assert names == {"cell_type", "disease"}
