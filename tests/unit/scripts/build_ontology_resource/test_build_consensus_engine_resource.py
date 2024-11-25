"""
Test for Cell Ontology Resource Module

This test suite validates the functionality of the `build_nx_graph_from_cl_ontology` and
`build_cell_ontology_ancestor_dictionary` functions in the Cell Ontology Resource module.

The tests cover various scenarios using both real and mocked `networkx` graphs to ensure the
correct processing of ontology data, including ancestor relationships.

Fixtures:
---------
- `cl_graph_meta_and_expected`: Generates real `networkx` graphs and expected ancestor dictionaries.

Test Cases:
-----------
- **Simple Chains**: Linear relationships (e.g., A -> B -> C).
- **Common Ancestors**: Multiple nodes sharing an ancestor.
- **Branching and Multiple Parents**: Handling of complex node relationships.
- **Edge Cases**: Empty graphs, isolated nodes, and nodes with no relationships.
"""

import typing as t

import networkx as nx
import pytest

from casp.scripts.build_ontology_resources.consensus_engine_resource import build_cell_ontology_ancestor_dictionary


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
    cl_graph_meta_and_expected: t.Tuple[nx.DiGraph, t.List[str], t.Dict[str, int], t.Dict[str, t.Any]]
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
