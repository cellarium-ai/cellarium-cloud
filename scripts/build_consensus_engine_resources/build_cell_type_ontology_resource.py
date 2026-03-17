import json
import logging
import typing as t

import networkx as nx
import owlready2
from smart_open import open

CL_PREFIX = "CL_"
logging.basicConfig(level=logging.INFO)


def build_nx_graph_from_cl_ontology(
    cl_ontology: owlready2.Ontology, cl_classes: t.List[owlready2.EntityClass]
) -> nx.DiGraph:
    """
    Build a directed graph from a cl ontology. Ontology object only considers immediate parents, while traversing the
    ontology ancestors, while in this algorithm we need to consider all ancestors. Network X library provides this
    functionality as we build the whole graph based on immediate parents and children of each node.

    :param cl_ontology: Cell Type Ontology object
    :param cl_classes: List of Cell Type Ontology classes

    :return: Directed graph
    """
    logging.info("Building nx graph from cl ontology...")
    cl_graph = nx.DiGraph(name="CL graph")
    cl_classes_set = set(cl_classes)

    for cl_class in cl_classes:
        cl_graph.add_node(cl_class.name)

    for cl_class in cl_classes:
        for parent_cl_class in cl_ontology.get_parents_of(entity=cl_class):
            if parent_cl_class in cl_classes_set:
                cl_graph.add_edge(u_of_edge=parent_cl_class.name, v_of_edge=cl_class.name)

        for child_cl_class in cl_ontology.get_children_of(Class=cl_class):
            if child_cl_class in cl_classes_set:
                cl_graph.add_edge(u_of_edge=cl_class.name, v_of_edge=child_cl_class.name)

    return cl_graph


def build_cell_ontology_ancestor_dictionary(
    cl_graph: nx.DiGraph, cl_names: t.List[str], cl_names_to_idx_map: t.Dict[str, int]
) -> t.Dict[str, t.List[str]]:
    """
    Build cell ontology ancestor dictionary resource. This resource allows a fast indexing of each CL node. It
    contains all ancestors for each of the CL node. It allows to fetch all ancestors for any node with a runtime O(1).

    :param cl_graph: CL ontology graph
    :param cl_names: list of CL names
    :param cl_names_to_idx_map: map from CL name to index
    """
    logging.info("Building cell ontology ancestor dictionary...")
    cl_ancestors_dictionary = {}

    for cl_name in cl_names:
        current_ancestors = []
        orders = []

        for cl_ancestor_name in sorted(nx.ancestors(cl_graph, cl_name)):
            cell_order_idx = cl_names_to_idx_map[cl_ancestor_name]
            current_ancestors.append(cl_ancestor_name)
            orders.append(cell_order_idx)

        sorted_ancestors = [ancestor for _, ancestor in sorted(zip(orders, current_ancestors))]

        cl_ancestors_dictionary[cl_name] = sorted_ancestors

    return cl_ancestors_dictionary


def main(cell_type_ontology_owl_file_url: str, output_file_path: str) -> None:
    """
    Create Cell Type Ontology Resources used in Ontology Aware Strategy and save it as a JSON file.

    :param cell_type_ontology_owl_file_url: Url for the OWL file to be used for the Cell Type Ontology
    :param output_file_path: Path to the output file (local or GCS)
    """
    logging.info("Generating cell type ontology resource from CL ontology...")
    cl_ontology = owlready2.get_ontology(cell_type_ontology_owl_file_url).load()

    cl_classes = [
        _class for _class in cl_ontology.classes() if _class.name.startswith(CL_PREFIX) and len(_class.label) == 1
    ]

    cl_names = list(set(_class.name for _class in cl_classes))
    cl_labels = list(set(_class.label[0] for _class in cl_classes))

    if len(cl_labels) != len(cl_names):
        raise ValueError("Number of unique cl labels doesn't correspond to number of unique cl names")

    cl_names_to_labels_map = {cl_name: cl_label for cl_name, cl_label in zip(cl_names, cl_labels)}
    cl_names_to_idx_map = {cl_name: idx for idx, cl_name in enumerate(cl_names)}

    cl_graph = build_nx_graph_from_cl_ontology(cl_ontology=cl_ontology, cl_classes=cl_classes)
    cl_ancestors_dictionary = build_cell_ontology_ancestor_dictionary(
        cl_graph=cl_graph, cl_names=cl_names, cl_names_to_idx_map=cl_names_to_idx_map
    )
    cell_ontology_resource = {
        "ancestors_dictionary": cl_ancestors_dictionary,
        "cell_ontology_term_id_to_cell_type": cl_names_to_labels_map,
    }

    logging.info(f"Writing output file to {output_file_path}")
    with open(output_file_path, "w") as output_file:
        json.dump(cell_ontology_resource, output_file)
