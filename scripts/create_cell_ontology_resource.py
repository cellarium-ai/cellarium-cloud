"""
Create the Cell Type Ontology resource JSON used by the Cellarium CAS backend.

The output JSON contains five precomputed fields consumed by the API:

  - ancestors_dictionary
  - cell_ontology_term_id_to_cell_type
  - children_dictionary
  - shortest_path_lengths_from_cell_root
  - longest_path_lengths_from_cell_root

Example::

    python scripts/create_cell_ontology_resource.py \\
        --owl-url https://github.com/obophenotype/cell-ontology/releases/download/v2025-07-30/cl.owl \\
        --output gs://your-bucket/path/to/cell_ontology_resource.json
"""

import json
import logging
import typing as t

import click
import networkx as nx
import owlready2
from smart_open import open

CELL_ROOT = "CL:0000000"
_OWL_NAME_PREFIX = "CL_"

logging.basicConfig(level=logging.INFO)


def _owl_name_to_cl_id(owl_name: str) -> str:
    """Convert owlready2 internal name (CL_XXXXXXX) to standard CL ID (CL:XXXXXXX)."""
    return owl_name.replace("CL_", "CL:", 1)


def build_nx_graph_from_cl_ontology(
    cl_ontology: owlready2.Ontology, cl_classes: t.List[owlready2.EntityClass]
) -> nx.DiGraph:
    """
    Build a directed graph from a CL ontology. Edges go from parent to child.
    Node names use standard CL: ID format.
    """
    logging.info("Building nx graph from cl ontology...")
    cl_graph = nx.DiGraph(name="CL graph")
    cl_classes_set = set(cl_classes)

    for cl_class in cl_classes:
        cl_graph.add_node(_owl_name_to_cl_id(cl_class.name))

    for cl_class in cl_classes:
        cl_id = _owl_name_to_cl_id(cl_class.name)
        for parent_cl_class in cl_ontology.get_parents_of(entity=cl_class):
            if parent_cl_class in cl_classes_set:
                cl_graph.add_edge(_owl_name_to_cl_id(parent_cl_class.name), cl_id)
        for child_cl_class in cl_ontology.get_children_of(Class=cl_class):
            if child_cl_class in cl_classes_set:
                cl_graph.add_edge(cl_id, _owl_name_to_cl_id(child_cl_class.name))

    return cl_graph


def build_ancestors_dictionary(
    cl_graph: nx.DiGraph, cl_names: t.List[str], cl_names_to_idx_map: t.Dict[str, int]
) -> t.Dict[str, t.List[str]]:
    """
    Build ancestor dictionary mapping each CL term to an ordered list of all its ancestors.
    Ordering follows the original index ordering of cl_names.
    """
    logging.info("Building cell ontology ancestor dictionary...")
    ancestors_dictionary = {}

    for cl_name in cl_names:
        current_ancestors = list(sorted(nx.ancestors(cl_graph, cl_name)))
        ancestors_dictionary[cl_name] = [
            a for _, a in sorted(zip([cl_names_to_idx_map[a] for a in current_ancestors], current_ancestors))
        ]

    return ancestors_dictionary


def build_children_dictionary(cl_graph: nx.DiGraph, cl_names: t.List[str]) -> t.Dict[str, t.List[str]]:
    """
    Build children dictionary mapping each CL term to its direct (non-transitive) children.
    """
    logging.info("Building cell ontology children dictionary...")
    return {cl_name: sorted(cl_graph.successors(cl_name)) for cl_name in cl_names}


def compute_shortest_path_lengths_from_root(cl_graph: nx.DiGraph) -> t.Dict[str, int]:
    """
    Compute the shortest path length from CL:0000000 to every reachable descendant. Root = 0.
    """
    logging.info("Computing shortest path lengths from root...")
    return dict(nx.single_source_shortest_path_length(cl_graph, CELL_ROOT))


def compute_longest_path_lengths_from_root(cl_graph: nx.DiGraph) -> t.Dict[str, int]:
    """
    Compute the longest path length from CL:0000000 to every reachable descendant using
    topological-order dynamic programming. Root = 0. Unreachable nodes are excluded.
    """
    logging.info("Computing longest path lengths from root...")
    dist: t.Dict[str, int] = {}

    for node in nx.topological_sort(cl_graph):
        if node == CELL_ROOT:
            dist[node] = 0
        if node not in dist:
            continue
        for child in cl_graph.successors(node):
            if child not in dist or dist[child] < dist[node] + 1:
                dist[child] = dist[node] + 1

    return dist


@click.command()
@click.option(
    "--owl-url",
    required=True,
    help="URL or local path to the CL OWL file.",
)
@click.option(
    "--output",
    required=True,
    help="Output file path (local or GCS, e.g. gs://bucket/path/resource.json).",
)
def main(owl_url: str, output: str) -> None:
    """Create the Cell Type Ontology resource JSON for Cellarium CAS."""
    logging.info("Generating cell type ontology resource from CL ontology...")

    # cl.owl imports ontologies (e.g. STATO) that define entities as both a property and
    # a class, which owlready2 rejects in _load_properties. Those entities are outside the
    # CL namespace and do not affect output, so we silence the specific TypeError.
    _orig = owlready2.Ontology._load_properties

    def _safe_load_properties(self: owlready2.Ontology) -> None:
        try:
            _orig(self)
        except TypeError as exc:
            logging.warning("Skipping property-loading error (imported ontology): %s", exc)

    owlready2.Ontology._load_properties = _safe_load_properties
    try:
        cl_ontology = owlready2.get_ontology(owl_url).load()
    finally:
        owlready2.Ontology._load_properties = _orig

    cl_classes = [
        _class
        for _class in cl_ontology.classes()
        if _class.name.startswith(_OWL_NAME_PREFIX) and len(_class.label) == 1
    ]

    # Build id->label mapping directly from class objects to avoid set-ordering bugs
    seen: t.Dict[str, str] = {}
    for _class in cl_classes:
        cl_id = _owl_name_to_cl_id(_class.name)
        if cl_id not in seen:
            seen[cl_id] = _class.label[0]

    cl_names_to_labels_map = seen
    cl_names = sorted(cl_names_to_labels_map.keys())
    cl_names_to_idx_map = {cl_name: idx for idx, cl_name in enumerate(cl_names)}

    cl_graph = build_nx_graph_from_cl_ontology(cl_ontology=cl_ontology, cl_classes=cl_classes)

    cell_ontology_resource = {
        "ancestors_dictionary": build_ancestors_dictionary(
            cl_graph=cl_graph, cl_names=cl_names, cl_names_to_idx_map=cl_names_to_idx_map
        ),
        "cell_ontology_term_id_to_cell_type": cl_names_to_labels_map,
        "children_dictionary": build_children_dictionary(cl_graph=cl_graph, cl_names=cl_names),
        "shortest_path_lengths_from_cell_root": compute_shortest_path_lengths_from_root(cl_graph=cl_graph),
        "longest_path_lengths_from_cell_root": compute_longest_path_lengths_from_root(cl_graph=cl_graph),
    }

    logging.info(f"Writing output file to {output}")
    with open(output, "w") as output_file:
        json.dump(cell_ontology_resource, output_file)


if __name__ == "__main__":
    main()
