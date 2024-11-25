"""
Test for Ontology and Graph Utilities

This test suite validates the functionality of various utility functions in the ontology processing module,
including the construction of NetworkX graphs from ontologies and the retrieval of ancestors and descendants
in directed graphs.

The tests use both mocked ontology data (via `owlready2` and `MagicMock`) and real `networkx` graphs to
ensure accurate representation and traversal of relationships within the ontology.

Fixtures:
---------
- `mock_ontology`: Creates a mocked ontology based on parameterized node and edge data. It generates various
  ontology structures such as simple chains, branching, multiple parents, and edge cases like empty graphs
  and isolated nodes.
- `simple_graph`: Provides a moderately complex directed graph structure to validate ancestor and descendant
  traversal functions in a non-mocked environment.

Test Cases:
-----------
- **Graph Construction from Ontology**: Verifies that a directed graph is correctly constructed from a mocked
  ontology using immediate parent-child relationships.
- **N-Level Ancestors**: Tests the ability to retrieve ancestors up to `n` levels in both simple and branching
  graphs.
- **N-Level Descendants**: Tests the ability to retrieve descendants up to `n` levels, including cases with
  branching nodes.
- **All Ancestors**: Validates that all ancestors of a given node, including transitive ancestors, are correctly
  retrieved.
- **All Descendants**: Ensures that all descendants of a given node, including transitive descendants, are
  correctly retrieved.
"""

from unittest.mock import MagicMock

import networkx as nx
import owlready2
import pytest

from casp.scripts.build_ontology_resources import utils


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


@pytest.fixture
def simple_graph() -> nx.DiGraph:
    """
    Create a simple directed graph for testing with moderate complexity.

    Graph structure:
        A -> B -> C
        A -> D
        D -> E
    """
    graph = nx.DiGraph()
    graph.add_edges_from([("A", "B"), ("B", "C"), ("A", "D"), ("D", "E")])
    return graph


def test_build_nx_graph_from_cl_ontology(mock_ontology: MagicMock) -> None:
    """
    Test the construction of a NetworkX graph from an ontology mock.

    This test verifies that the graph structure and relationships are accurately represented when building the graph
    using the mocked ontology. It checks that nodes, ancestors, and descendants are correctly identified.

    :param mock_ontology: The mocked ontology object created by the fixture.
    """
    cl_classes = list(mock_ontology.classes())

    # Build the graph using the mocked ontology with the real method replacements
    cl_graph = utils.build_nx_graph_from_cl_ontology(cl_ontology=mock_ontology, cl_classes=cl_classes)

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


def test_get_n_level_ancestors(simple_graph: nx.DiGraph) -> None:
    """
    Test get_n_level_ancestors function for various nodes and levels.
    """
    assert utils.get_n_level_ancestors(simple_graph, "C", 1) == {"B", "C"}, "Incorrect 1-level ancestors for node 'C'."
    assert utils.get_n_level_ancestors(simple_graph, "C", 2) == {
        "A",
        "B",
        "C",
    }, "Incorrect 2-level ancestors for node 'C'."
    assert utils.get_n_level_ancestors(simple_graph, "E", 1) == {"D", "E"}, "Incorrect 1-level ancestors for node 'E'."
    assert utils.get_n_level_ancestors(simple_graph, "E", 2) == {
        "A",
        "D",
        "E",
    }, "Incorrect 2-level ancestors for node 'E'."


def test_get_n_level_descendants(simple_graph: nx.DiGraph) -> None:
    """
    Test get_n_level_descendants function for various nodes and levels.
    """
    assert utils.get_n_level_descendants(simple_graph, "A", 1) == {
        "A",
        "B",
        "D",
    }, "Incorrect 1-level descendants for node 'A'."
    assert utils.get_n_level_descendants(simple_graph, "A", 2) == {
        "A",
        "B",
        "C",
        "D",
        "E",
    }, "Incorrect 2-level descendants for node 'A'."
    assert utils.get_n_level_descendants(simple_graph, "B", 1) == {
        "B",
        "C",
    }, "Incorrect 1-level descendants for node 'B'."
    assert utils.get_n_level_descendants(simple_graph, "D", 1) == {
        "D",
        "E",
    }, "Incorrect 1-level descendants for node 'D'."


def test_get_all_ancestors(simple_graph: nx.DiGraph) -> None:
    """
    Test get_all_ancestors function for various nodes.
    """
    assert utils.get_all_ancestors(simple_graph, "C") == {"A", "B", "C"}, "Incorrect ancestors for node 'C'."
    assert utils.get_all_ancestors(simple_graph, "E") == {"A", "D", "E"}, "Incorrect ancestors for node 'E'."
    assert utils.get_all_ancestors(simple_graph, "A") == {"A"}, "Incorrect ancestors for node 'A'."


def test_get_all_descendants(simple_graph: nx.DiGraph) -> None:
    """
    Test get_all_descendants function for various nodes.
    """
    assert utils.get_all_descendants(simple_graph, "A") == {
        "A",
        "B",
        "C",
        "D",
        "E",
    }, "Incorrect descendants for node 'A'."
    assert utils.get_all_descendants(simple_graph, "B") == {"B", "C"}, "Incorrect descendants for node 'B'."
    assert utils.get_all_descendants(simple_graph, "D") == {"D", "E"}, "Incorrect descendants for node 'D'."
