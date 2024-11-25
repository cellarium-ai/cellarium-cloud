import json
import logging
import typing as t

import networkx as nx
import owlready2
from smart_open import open

from casp.scripts.build_ontology_resources import utils

CL_PREFIX = "CL_"
logging.basicConfig(level=logging.INFO)


def build_benchmarking_ontology_dictionary_resource(
    cl_graph: nx.DiGraph, cl_names: t.List[str], n_hops: int
) -> t.Dict[str, t.Any]:
    """
    Build a cell ontology ancestor and descendant dictionary resource.

    This resource provides fast indexing of each CL node up to a specified number of hops (`n_hops`).
    It includes ancestors and descendants for all nodes, along with hop level information.
    The hop level information contains:

        * Ancestor nodes for each hop, up to `n_hops` levels deep:
            - If `hop_1`, it includes immediate ancestors.
            - If `hop_2`, it includes ancestors of ancestors, and so on.

        * Descendant nodes for each hop, up to `n_hops` levels deep:
            - If `hop_1`, it includes immediate descendants.
            - If `hop_2`, it includes descendants of descendants, and so on.

        * Union of all ancestor and descendant nodes for each hop level.

    The resource allows retrieval of all ancestors and descendants for any node at any hop with O(1) runtime.

    Output Example
    --------------

    .. code-block:: python

        {
            "key1": "value1",
            "key2": "value2",
            "nested_key": {
                "subkey1": "subvalue1",
                "subkey2": "subvalue2"
            }
        }

    :param cl_graph: The CL ontology graph.
    :param cl_names: A list of CL names.
    :param n_hops: The number of hops to consider.
    """
    ontology_resource_dict = {}

    for cl_name in cl_names:
        ontology_resource_dict[cl_name] = {}
        ontology_resource_dict[cl_name]["all_ancestors"] = utils.get_all_ancestors(cl_graph, node=cl_name)
        ontology_resource_dict[cl_name]["all_descendants"] = utils.get_all_descendants(cl_graph, node=cl_name)

        for top_n in range(n_hops):
            nodes = utils.get_n_level_ancestors(cl_graph, node=cl_name, n=top_n)
            hop_all_ancestors = set.union(
                *[utils.get_all_ancestors(cl_graph, node=node_cl_name) for node_cl_name in list(nodes)]
            )
            hop_all_descendants = set.union(
                *[utils.get_all_descendants(cl_graph, node=node_cl_name) for node_cl_name in list(nodes)]
            )

            ontology_resource_dict[cl_name][f"hop_{top_n}"] = {
                "nodes": nodes,
                "all_ancestors": hop_all_ancestors,
                "all_descendants": hop_all_descendants,
            }

    return ontology_resource_dict


def main(cell_type_ontology_owl_file_url: str, output_file_path: str, n_hops: int) -> None:
    """
    Create Cell Type Ontology Resources used in Ontology Aware Strategy and save it as a JSON file.

    :param cell_type_ontology_owl_file_url: Url for the OWL file to be used for the Cell Type Ontology.
    :param output_file_path: Path to the output file (local or GCS).
    :param n_hops: Number of hops to include in the resource file.
    """
    logging.info("Generating cell type ontology resource for benchmarking from CL ontology...")
    cl_ontology = owlready2.get_ontology(cell_type_ontology_owl_file_url).load()

    cl_classes = [
        _class for _class in cl_ontology.classes() if _class.name.startswith(CL_PREFIX) and len(_class.label) == 1
    ]

    cl_names = list(set(_class.name for _class in cl_classes))
    cl_labels = list(set(_class.label[0] for _class in cl_classes))

    if len(cl_labels) != len(cl_names):
        raise ValueError("Number of unique cl labels doesn't correspond to number of unique cl names")

    cl_graph = utils.build_nx_graph_from_cl_ontology(cl_ontology=cl_ontology, cl_classes=cl_classes)

    ontology_resource_dict = build_benchmarking_ontology_dictionary_resource(
        cl_graph=cl_graph, cl_names=cl_names, n_hops=n_hops
    )

    logging.info(f"Writing output file to {output_file_path}")

    with open(output_file_path, "w") as output_file:
        json.dump(ontology_resource_dict, output_file)
