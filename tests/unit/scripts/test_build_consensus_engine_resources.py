"""
Test for Cell Ontology Resource Module

This test suite validates the functionality of the `build_nx_graph_from_cl_ontology` and
`build_cell_ontology_ancestor_dictionary` functions in the Cell Ontology Resource module.

The tests cover various scenarios using both real and mocked `networkx` graphs to ensure the
correct processing of ontology data, including ancestor relationships.

Fixtures:
---------
- `mock_ontology`: Creates mocked owlready2.Ontology
- `cl_graph_meta_and_expected`: Generates real `networkx` graphs and expected ancestor dictionaries.

Test Cases:
-----------
- **Simple Chains**: Linear relationships (e.g., A -> B -> C).
- **Common Ancestors**: Multiple nodes sharing an ancestor.
- **Branching and Multiple Parents**: Handling of complex node relationships.
- **Edge Cases**: Empty graphs, isolated nodes, and nodes with no relationships.
"""

import typing as t
from unittest.mock import MagicMock

import networkx as nx
import owlready2
import pytest

from casp.scripts.build_consensus_engine_resources.build_cell_type_ontology_resource import (
    build_cell_ontology_ancestor_dictionary,
    build_nx_graph_from_cl_ontology,
)


@pytest.fixture(
    params=[
        # Test case 1: Simple chain A -> B -> C
        {"nodes": ["A", "B", "C"], "edges": {"A": [], "B": ["A"], "C": ["B"]}},
        # Test case 2: Common ancestor A -> B, A -> C
        {"nodes": ["A", "B", "C"], "edges": {"A": [], "B": ["A"], "C": ["A"]}},
        # Test case 3: Branching A -> B -> C, A -> D
        {"nodes": ["A", "B", "C", "D"], "edges": {"A": [], "B": ["A"], "C": ["B"], "D": ["A"]}},
        # Test case 4: Multiple parents A -> B, B -> D, A -> C, C -> D
        {"nodes": ["A", "B", "C", "D"], "edges": {"A": [], "B": ["A"], "C": ["B"], "D": ["C", "B"]}},
        # Edge case 1: Empty graph (no nodes)
        {"nodes": [], "edges": {}},
        # Edge case 2: Single node with no parents or children
        {"nodes": ["A"], "edges": {"A": []}},
        # Edge case 3: Node with no ancestors or descendants
        # A -> B, C is isolated
        {"nodes": ["A", "B", "C"], "edges": {"A": [], "B": ["A"], "C": []}},
    ]
)
def mock_ontology(request: pytest.FixtureRequest) -> MagicMock:
    """
    Create a mocked ontology based on parameterized node and edge data.

    The fixture generates a mock ontology where nodes represent classes and edges represent parent-child relationships.
    It covers various scenarios, including simple chains, branching, multiple parents, and edge cases like empty graphs
    and isolated nodes.

    :param request: The pytest request object that provides access to the parameterized data.

    :return: A mocked ontology object with classes and relationships.
    """
    # Create mock classes based on the parameterized relationships
    class_map = {}
    for class_name in request.param["nodes"]:
        mock_class = MagicMock(spec=owlready2.EntityClass)
        mock_class.name = class_name  # Use the original class name as it is
        mock_class.label = [class_name]  # Use the class name directly as the label
        class_map[class_name] = mock_class

    # Create the mock ontology
    cl_ontology = MagicMock(spec=owlready2.Ontology)
    cl_ontology.classes.return_value = list(class_map.values())

    # Define the function to mock get_parents_of
    def get_parents_of(entity):
        return [class_map[parent] for parent in request.param["edges"].get(entity.name, [])]

    # Define the function to mock get_children_of
    def get_children_of(Class):
        return [class_map[child] for child, parents in request.param["edges"].items() if Class.name in parents]

    # Mock the get_parents_of and get_children_of methods
    cl_ontology.get_parents_of.side_effect = get_parents_of
    cl_ontology.get_children_of.side_effect = get_children_of

    return cl_ontology


def test_build_nx_graph_from_cl_ontology(mock_ontology: MagicMock) -> None:
    """
    Test the construction of a NetworkX graph from an ontology mock.

    This test verifies that the graph structure and relationships are accurately represented when building the graph
    using the mocked ontology. It checks that nodes, ancestors, and descendants are correctly identified.

    :param mock_ontology: The mocked ontology object created by the fixture.
    """
    cl_classes = list(mock_ontology.classes())

    # Build the graph using the mocked ontology with the real method replacements
    cl_graph = build_nx_graph_from_cl_ontology(cl_ontology=mock_ontology, cl_classes=cl_classes)

    # Verify the structure of the graph
    class_names = [cl.name for cl in cl_classes]
    assert set(cl_graph.nodes()) == set(class_names), "The nodes in the graph do not match the expected class names."

    # Check ancestors and descendants for each class
    for cl_class in cl_classes:
        name = cl_class.name
        parents = mock_ontology.get_parents_of(cl_class)
        ancestors = nx.ancestors(cl_graph, name)
        for parent in parents:
            assert parent.name in ancestors, f"{parent.name} is not present in ancestors of cl graph"

        children = mock_ontology.get_children_of(cl_class)
        descendants = nx.descendants(cl_graph, name)
        for child in children:
            assert child.name in descendants, f"{child.name} is not present in descendants of cl graph"


@pytest.fixture(
    params=[
        # Test case 1: Simple chain A -> B -> C
        {
            "nodes": ["A", "B", "C"],
            "edges": [("A", "B"), ("B", "C")],
            "expected": {"A": [], "B": ["A"], "C": ["A", "B"]},
        },
        # Test case 2: Common ancestor A -> B, A -> C
        {"nodes": ["A", "B", "C"], "edges": [("A", "B"), ("A", "C")], "expected": {"A": [], "B": ["A"], "C": ["A"]}},
        # Test case 3: Branching A -> B -> C, A -> D
        {
            "nodes": ["A", "B", "C", "D"],
            "edges": [("A", "B"), ("B", "C"), ("A", "D")],
            "expected": {"A": [], "B": ["A"], "C": ["A", "B"], "D": ["A"]},
        },
        # Test case 4: Multiple parents A -> B, B -> D, A -> C, C -> D
        {
            "nodes": ["A", "B", "C", "D"],
            "edges": [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")],
            "expected": {"A": [], "B": ["A"], "C": ["A"], "D": ["A", "B", "C"]},
        },
        # Edge case 1: Empty graph (no nodes)
        {"nodes": [], "edges": [], "expected": {}},
        # Edge case 2: Single node with no parents or children
        {"nodes": ["A"], "edges": [], "expected": {"A": []}},
        # Edge case 3: Node with no ancestors or descendants
        # A -> B, C is isolated
        {
            "nodes": ["A", "B", "C"],
            "edges": [("A", "B")],
            "expected": {"A": [], "B": ["A"], "C": []},  # C is isolated and has no ancestors
        },
    ]
)
def cl_graph_meta_and_expected(
    request: pytest.FixtureRequest,
) -> t.Tuple[nx.DiGraph, t.List[str], t.Dict[str, int], t.Dict[str, t.Any]]:
    """
    Fixture that generates a NetworkX graph and the expected ancestor dictionary based on parameterized node and edge
    data.

    This fixture creates a directed graph from the provided nodes and edges, ensuring all nodes are present, even if
    isolated. It also provides the expected output for the ancestor dictionary based on the graph structure.

    :param request: The pytest request object that provides access to the parameterized data.
    :return: A tuple containing the graph, list of class names, index map, and the expected ancestor dictionary.
    """
    # Create a directed graph and add edges
    cl_graph = nx.DiGraph()
    cl_graph.add_nodes_from(request.param["nodes"])
    cl_graph.add_edges_from(request.param["edges"])

    # Define class names and index map
    cl_names = list(cl_graph.nodes())
    cl_names_to_idx_map = {name: i for i, name in enumerate(cl_names)}

    return cl_graph, cl_names, cl_names_to_idx_map, request.param["expected"]


def test_build_cell_ontology_ancestor_dictionary(
    cl_graph_meta_and_expected: t.Tuple[nx.DiGraph, t.List[str], t.Dict[str, int], t.Dict[str, t.Any]],
):
    """
    Test the construction of the cell ontology ancestor dictionary from a NetworkX graph.

    This test verifies that the ancestor dictionary is correctly built from the graph, ensuring that all ancestors are
    accurately identified and ordered according to the index map. It compares the generated dictionary against the
    expected output.

    :param cl_graph_meta_and_expected: The tuple returned by the fixture, containing the graph, class names, index map,
        and expected output.
    """
    cl_graph, cl_names, cl_names_to_idx_map, expected = cl_graph_meta_and_expected

    # Call the function
    ancestor_dict = build_cell_ontology_ancestor_dictionary(cl_graph, cl_names, cl_names_to_idx_map)

    # Assert that the output matches the expected dictionary
    assert ancestor_dict == expected, "Output of `build_cell_ontology_ancestor_dictionary` is note as expected"
