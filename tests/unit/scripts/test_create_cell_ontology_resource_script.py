"""
Unit tests for scripts/create_cell_ontology_resource.py.

These tests exercise every function in the script using a small hand-crafted
networkx graph so that owlready2 and a real OWL file are not required.
"""

import json
from pathlib import Path
import tempfile
from unittest.mock import MagicMock, patch

from click.testing import CliRunner
import networkx as nx
import pytest

from scripts.create_cell_ontology_resource import (
    CELL_ROOT,
    _owl_name_to_cl_id,
    build_ancestors_dictionary,
    build_children_dictionary,
    build_nx_graph_from_cl_ontology,
    compute_longest_path_lengths_from_root,
    compute_shortest_path_lengths_from_root,
    main,
)


def _make_simple_graph() -> nx.DiGraph:
    """
    Simple 4-node DAG:
        CL:0000000 → CL:0000001 → CL:0000002
                   ↘              ↗
                    CL:0000003
    CL:0000002 is reachable from root via two paths of lengths 2 and 2,
    but the longest is 2 via CL:0000003 too — both paths len 2.
    """
    g = nx.DiGraph()
    g.add_edges_from(
        [
            ("CL:0000000", "CL:0000001"),
            ("CL:0000001", "CL:0000002"),
            ("CL:0000000", "CL:0000003"),
            ("CL:0000003", "CL:0000002"),
        ]
    )
    return g


def _make_deep_graph() -> nx.DiGraph:
    """Linear chain: root → A → B → C"""
    g = nx.DiGraph()
    g.add_edges_from(
        [
            (CELL_ROOT, "CL:0000001"),
            ("CL:0000001", "CL:0000002"),
            ("CL:0000002", "CL:0000003"),
        ]
    )
    return g


def test_owl_name_to_cl_id_replaces_prefix():
    """Standard CL_ prefix is converted to CL:."""
    assert _owl_name_to_cl_id("CL_0000001") == "CL:0000001"


def test_owl_name_to_cl_id_only_replaces_first_occurrence():
    """Only the leading CL_ prefix is replaced; subsequent occurrences are left intact."""
    assert _owl_name_to_cl_id("CL_CL_0000001") == "CL:CL_0000001"


def test_owl_name_to_cl_id_non_cl_name_unchanged():
    """Non-CL ontology names are returned without modification."""
    assert _owl_name_to_cl_id("GO_0008150") == "GO_0008150"


def _make_mock_class(name: str, label: str, parents=None, children=None):
    cls = MagicMock()
    cls.name = name
    cls.label = [label]
    cls.__hash__ = lambda self: hash(name)
    cls.__eq__ = lambda self, other: self.name == getattr(other, "name", None)
    return cls


def test_build_nx_graph_nodes_use_cl_colon_format():
    """Graph nodes must use CL: format, not the owlready2 internal CL_ format."""
    a = _make_mock_class("CL_0000000", "cell")
    b = _make_mock_class("CL_0000001", "somatic cell")

    ontology = MagicMock()
    ontology.get_parents_of.side_effect = lambda entity: [a] if entity is b else []
    ontology.get_children_of.side_effect = lambda **kwargs: [b] if kwargs["Class"] is a else []

    graph = build_nx_graph_from_cl_ontology(cl_ontology=ontology, cl_classes=[a, b])

    assert "CL:0000000" in graph.nodes
    assert "CL:0000001" in graph.nodes


def test_build_nx_graph_edges_parent_to_child():
    """Edges are directed parent → child, not child → parent."""
    a = _make_mock_class("CL_0000000", "cell")
    b = _make_mock_class("CL_0000001", "somatic cell")

    ontology = MagicMock()
    ontology.get_parents_of.side_effect = lambda entity: [a] if entity is b else []
    ontology.get_children_of.side_effect = lambda **kwargs: [b] if kwargs["Class"] is a else []

    graph = build_nx_graph_from_cl_ontology(cl_ontology=ontology, cl_classes=[a, b])

    assert graph.has_edge("CL:0000000", "CL:0000001")
    assert not graph.has_edge("CL:0000001", "CL:0000000")


def test_build_nx_graph_excludes_out_of_set_parents():
    """Classes outside the filtered CL set are ignored and create no edges."""
    a = _make_mock_class("CL_0000000", "cell")
    b = _make_mock_class("CL_0000001", "somatic cell")
    external = _make_mock_class("GO_0000000", "external")

    ontology = MagicMock()
    ontology.get_parents_of.side_effect = lambda entity: [external] if entity is b else []
    ontology.get_children_of.return_value = []

    graph = build_nx_graph_from_cl_ontology(cl_ontology=ontology, cl_classes=[a, b])

    assert graph.number_of_edges() == 0


def test_build_ancestors_dictionary_root_has_no_ancestors():
    """The ontology root (CL:0000000) has no ancestors."""
    g = _make_simple_graph()
    cl_names = sorted(g.nodes)
    idx_map = {n: i for i, n in enumerate(cl_names)}
    result = build_ancestors_dictionary(g, cl_names, idx_map)
    assert result["CL:0000000"] == []


def test_build_ancestors_dictionary_leaf_includes_all_ancestors():
    """Every node on the path from root to a leaf appears in its ancestor list."""
    g = _make_deep_graph()
    cl_names = sorted(g.nodes)
    idx_map = {n: i for i, n in enumerate(cl_names)}
    result = build_ancestors_dictionary(g, cl_names, idx_map)
    assert set(result["CL:0000003"]) == {"CL:0000000", "CL:0000001", "CL:0000002"}


def test_build_ancestors_dictionary_ordering_follows_idx_map():
    """Ancestors are returned in ascending cl_names index order, not arbitrary set order."""
    g = _make_deep_graph()
    cl_names = sorted(g.nodes)
    idx_map = {n: i for i, n in enumerate(cl_names)}
    result = build_ancestors_dictionary(g, cl_names, idx_map)
    ancestors = result["CL:0000003"]
    assert ancestors == [cl_names[idx_map[a]] for a in sorted(ancestors, key=lambda x: idx_map[x])]


def test_build_ancestors_dictionary_all_keys_present():
    """Every CL term in the graph has an entry in the output dictionary."""
    g = _make_simple_graph()
    cl_names = sorted(g.nodes)
    idx_map = {n: i for i, n in enumerate(cl_names)}
    result = build_ancestors_dictionary(g, cl_names, idx_map)
    assert set(result.keys()) == set(cl_names)


def test_build_children_dictionary_root_has_correct_children():
    """Root's direct children are the two nodes immediately below it in the DAG."""
    g = _make_simple_graph()
    cl_names = sorted(g.nodes)
    result = build_children_dictionary(g, cl_names)
    assert set(result["CL:0000000"]) == {"CL:0000001", "CL:0000003"}


def test_build_children_dictionary_leaf_has_empty_list():
    """Leaf nodes have an empty children list, not a missing key."""
    g = _make_simple_graph()
    cl_names = sorted(g.nodes)
    result = build_children_dictionary(g, cl_names)
    assert result["CL:0000002"] == []


def test_build_children_dictionary_values_are_sorted():
    """Children lists are lexicographically sorted for deterministic output."""
    g = _make_simple_graph()
    cl_names = sorted(g.nodes)
    result = build_children_dictionary(g, cl_names)
    for children in result.values():
        assert children == sorted(children)


def test_build_children_dictionary_all_keys_present():
    """Every CL term in the graph has an entry, including internal and leaf nodes."""
    g = _make_deep_graph()
    cl_names = sorted(g.nodes)
    result = build_children_dictionary(g, cl_names)
    assert set(result.keys()) == set(cl_names)


def test_shortest_path_root_is_zero():
    """CL:0000000 has a path length of 0 from itself."""
    g = _make_deep_graph()
    result = compute_shortest_path_lengths_from_root(g)
    assert result[CELL_ROOT] == 0


def test_shortest_path_linear_chain():
    """Each hop along a linear chain increments the path length by 1."""
    g = _make_deep_graph()
    result = compute_shortest_path_lengths_from_root(g)
    assert result["CL:0000001"] == 1
    assert result["CL:0000002"] == 2
    assert result["CL:0000003"] == 3


def test_shortest_path_diamond_takes_shorter():
    """
    root → A (len 1)
    root → B → A is shorter via root→A directly.
    Diamond: root→A→C and root→B→C — shortest to C is 2.
    """
    g = nx.DiGraph()
    g.add_edges_from(
        [
            (CELL_ROOT, "CL:0000001"),  # root → A, len 1
            (CELL_ROOT, "CL:0000002"),  # root → B, len 1
            ("CL:0000001", "CL:0000003"),  # A → C, len 2
            ("CL:0000002", "CL:0000003"),  # B → C, len 2
        ]
    )
    result = compute_shortest_path_lengths_from_root(g)
    assert result["CL:0000003"] == 2


def test_shortest_path_excludes_unreachable_nodes():
    """Nodes with no path from the root are absent from the result."""
    g = _make_deep_graph()
    g.add_node("CL:9999999")  # disconnected
    result = compute_shortest_path_lengths_from_root(g)
    assert "CL:9999999" not in result


def test_longest_path_root_is_zero():
    """CL:0000000 has a longest path length of 0 from itself."""
    g = _make_deep_graph()
    result = compute_longest_path_lengths_from_root(g)
    assert result[CELL_ROOT] == 0


def test_longest_path_linear_chain():
    """On a linear chain, longest path equals shortest path at each node."""
    g = _make_deep_graph()
    result = compute_longest_path_lengths_from_root(g)
    assert result["CL:0000003"] == 3


def test_longest_path_prefers_longer_route():
    """
    root → A (len 1) → C (len 2)
    root → B (len 1) → A (len 2) → C (len 3)
    Longest path to C should be 3.
    """
    g = nx.DiGraph()
    g.add_edges_from(
        [
            (CELL_ROOT, "CL:0000001"),  # root → A
            (CELL_ROOT, "CL:0000002"),  # root → B
            ("CL:0000002", "CL:0000001"),  # B → A (makes A reachable at depth 2)
            ("CL:0000001", "CL:0000003"),  # A → C
        ]
    )
    result = compute_longest_path_lengths_from_root(g)
    assert result["CL:0000003"] == 3


def test_longest_path_excludes_unreachable_nodes():
    """Nodes with no path from the root are absent from the result."""
    g = _make_deep_graph()
    g.add_node("CL:9999999")  # disconnected
    result = compute_longest_path_lengths_from_root(g)
    assert "CL:9999999" not in result


def _make_mock_cl_class(name: str, label: str) -> MagicMock:
    obj = MagicMock()
    obj.name = name
    obj.label = [label]
    obj.__hash__ = lambda self: hash(name)
    obj.__eq__ = lambda self, other: getattr(other, "name", None) == name
    return obj


@pytest.fixture
def minimal_ontology_mock():
    """
    Three-node ontology: root → somatic → T cell (linear chain).
    CL_0000000 → CL_0000003 → CL_0000084
    """
    root = _make_mock_cl_class("CL_0000000", "cell")
    somatic = _make_mock_cl_class("CL_0000003", "somatic cell")
    tcell = _make_mock_cl_class("CL_0000084", "T cell")
    classes = [root, somatic, tcell]

    def get_parents_of(entity):
        if entity is somatic:
            return [root]
        if entity is tcell:
            return [somatic]
        return []

    def get_children_of(**kwargs):
        cl = kwargs["Class"]
        if cl is root:
            return [somatic]
        if cl is somatic:
            return [tcell]
        return []

    ontology = MagicMock()
    ontology.classes.return_value = iter(classes)
    ontology.get_parents_of.side_effect = get_parents_of
    ontology.get_children_of.side_effect = get_children_of
    return ontology


def test_main_writes_valid_json(minimal_ontology_mock):
    """CLI produces a JSON file containing exactly the five expected top-level keys."""
    runner = CliRunner()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = str(Path(tmpdir) / "resource.json")
        with patch("scripts.create_cell_ontology_resource.owlready2") as mock_owlready2:
            mock_owlready2.get_ontology.return_value.load.return_value = minimal_ontology_mock
            result = runner.invoke(main, ["--owl-url", "fake://cl.owl", "--output", output_path])

        assert result.exit_code == 0, result.output
        with open(output_path) as f:
            data = json.load(f)

    assert set(data.keys()) == {
        "ancestors_dictionary",
        "cell_ontology_term_id_to_cell_type",
        "children_dictionary",
        "shortest_path_lengths_from_cell_root",
        "longest_path_lengths_from_cell_root",
    }


def test_main_output_uses_cl_colon_format(minimal_ontology_mock):
    """All ontology term keys in the output use CL: format, not the owlready2 CL_ format."""
    runner = CliRunner()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = str(Path(tmpdir) / "resource.json")
        with patch("scripts.create_cell_ontology_resource.owlready2") as mock_owlready2:
            mock_owlready2.get_ontology.return_value.load.return_value = minimal_ontology_mock
            runner.invoke(main, ["--owl-url", "fake://cl.owl", "--output", output_path])

        with open(output_path) as f:
            data = json.load(f)

    for key in data["ancestors_dictionary"]:
        assert key.startswith("CL:"), f"Found non-CL: key: {key}"
    for key in data["cell_ontology_term_id_to_cell_type"]:
        assert key.startswith("CL:"), f"Found non-CL: key: {key}"


def test_main_ancestors_of_root_are_empty(minimal_ontology_mock):
    """The ontology root has no ancestors in the output JSON."""
    runner = CliRunner()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = str(Path(tmpdir) / "resource.json")
        with patch("scripts.create_cell_ontology_resource.owlready2") as mock_owlready2:
            mock_owlready2.get_ontology.return_value.load.return_value = minimal_ontology_mock
            runner.invoke(main, ["--owl-url", "fake://cl.owl", "--output", output_path])

        with open(output_path) as f:
            data = json.load(f)

    assert data["ancestors_dictionary"]["CL:0000000"] == []


def test_main_shortest_and_longest_agree_on_linear_chain(minimal_ontology_mock):
    """On a linear chain (no branching), shortest and longest path lengths are identical."""
    runner = CliRunner()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = str(Path(tmpdir) / "resource.json")
        with patch("scripts.create_cell_ontology_resource.owlready2") as mock_owlready2:
            mock_owlready2.get_ontology.return_value.load.return_value = minimal_ontology_mock
            runner.invoke(main, ["--owl-url", "fake://cl.owl", "--output", output_path])

        with open(output_path) as f:
            data = json.load(f)

    # Linear chain — shortest == longest for every node
    shortest = data["shortest_path_lengths_from_cell_root"]
    longest = data["longest_path_lengths_from_cell_root"]
    for key in shortest:
        assert shortest[key] == longest[key], f"Mismatch for {key}"


def test_main_missing_required_args_fails():
    """Invoking CLI without required arguments exits with a non-zero status."""
    runner = CliRunner()
    result = runner.invoke(main, [])
    assert result.exit_code != 0
