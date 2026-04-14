"""
Tests things in the cellarium_general_service that can reasonably be tested with a unit test.
"""

import json
import typing as t

from mockito import mock, unstub
import pytest

from cellarium.cas_backend.apps.compute.schemas import CellOntologyResourceResponse
from cellarium.cas_backend.apps.compute.services.cellarium_general_service import CellariumGeneralService
from cellarium.cas_backend.apps.compute.services.consensus_engine.strategies.ontology_aware import CellOntologyResource
from cellarium.cas_backend.apps.compute.services.exceptions import InvalidClientVersionException
from cellarium.cas_backend.core import constants


class TestCellariumGeneralService:
    def setup_method(self) -> None:
        self.cellarium_general_service = CellariumGeneralService(cellarium_general_dm=mock())

    def teardown_method(self) -> None:
        unstub()

    def test_validate_client_version_valid(self):
        assert self.cellarium_general_service.validate_client_version(constants.MIN_CLIENT_VERSION)

    def test_validate_client_version_too_old(self):
        assert not self.cellarium_general_service.validate_client_version("1.3.0")

    def test_validate_client_version_invalid(self):
        with pytest.raises(InvalidClientVersionException):
            self.cellarium_general_service.validate_client_version("invalid_version")


def _load_ontology_resource_dict() -> dict[str, t.Any]:
    filepath = "tests/unit/test_consensus_engine_fixtures/cell_ontology_resource_mini.json"
    with open(filepath) as f:
        return json.load(f)


class TestCellOntologyResource:
    def setup_method(self) -> None:
        self.resource_dict = _load_ontology_resource_dict()
        self.resource = CellOntologyResource(cell_ontology_resource_dict=self.resource_dict)

    def test_cl_names_is_sorted(self):
        cl_names = self.resource.cl_names
        assert cl_names == sorted(cl_names)

    def test_cl_names_contains_all_term_ids(self):
        cl_names = self.resource.cl_names
        assert set(cl_names) == set(self.resource_dict["cell_ontology_term_id_to_cell_type"].keys())

    def test_children_dictionary_loaded(self):
        assert self.resource.children_dictionary is not None
        assert self.resource.children_dictionary == self.resource_dict["children_dictionary"]

    def test_shortest_path_lengths_loaded(self):
        assert self.resource.shortest_path_lengths_from_cell_root is not None
        assert self.resource.shortest_path_lengths_from_cell_root["CL_0000000"] == 0

    def test_longest_path_lengths_loaded(self):
        assert self.resource.longest_path_lengths_from_cell_root is not None
        assert self.resource.longest_path_lengths_from_cell_root["CL_0000000"] == 0

    def test_new_fields_optional_when_missing(self):
        minimal_dict = {
            "ancestors_dictionary": self.resource_dict["ancestors_dictionary"],
            "cell_ontology_term_id_to_cell_type": self.resource_dict["cell_ontology_term_id_to_cell_type"],
        }
        resource = CellOntologyResource(cell_ontology_resource_dict=minimal_dict)
        assert resource.children_dictionary is None
        assert resource.shortest_path_lengths_from_cell_root is None
        assert resource.longest_path_lengths_from_cell_root is None


class TestGetCellOntologyResourceEndpoint:
    def setup_method(self) -> None:
        self.resource_dict = _load_ontology_resource_dict()

    def test_response_has_correct_keys(self):
        resource = CellOntologyResource(cell_ontology_resource_dict=self.resource_dict)
        response = CellOntologyResourceResponse(
            cl_names=resource.cl_names,
            cell_ontology_term_id_to_cell_type=resource.ontology_term_id_to_name_dict,
            children_dictionary=resource.children_dictionary,
            shortest_path_lengths_from_cell_root=resource.shortest_path_lengths_from_cell_root,
            longest_path_lengths_from_cell_root=resource.longest_path_lengths_from_cell_root,
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

    def test_cl_names_sorted_in_response(self):
        resource = CellOntologyResource(cell_ontology_resource_dict=self.resource_dict)
        response = CellOntologyResourceResponse(
            cl_names=resource.cl_names,
            cell_ontology_term_id_to_cell_type=resource.ontology_term_id_to_name_dict,
            children_dictionary=resource.children_dictionary,
            shortest_path_lengths_from_cell_root=resource.shortest_path_lengths_from_cell_root,
            longest_path_lengths_from_cell_root=resource.longest_path_lengths_from_cell_root,
        )
        assert response.cl_names == sorted(response.cl_names)
