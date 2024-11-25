import logging
import typing as t

import networkx as nx
import owlready2


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


def get_n_level_ancestors(graph: nx.DiGraph, node: str, n: int) -> t.Set[str]:
    """
    Return the set of all ancestors of a given node in a directed graph up to `n` levels deep.

    This function traverses `n` levels of the graph in reverse (from child to parent), starting from the specified
    `node`. It collects all nodes that are reachable from `node` within `n` levels of predecessors (ancestors).
    The original node is also included in the resulting set.

    :param graph: A directed graph represented by a NetworkX DiGraph object, where nodes represent entities and
        directed edges represent relationships between them.
    :param node: The node for which to find ancestors. This is the starting point of the traversal.
    :param n: The number of levels to traverse upwards in the graph, i.e., how far back to go in the graph to find
        ancestors. If n = 1, it returns immediate predecessors; if n = 2, it includes predecessors of predecessors,
        and so on.

    :return: A set containing the node itself and all ancestors reachable within `n` levels from the node. If the node
        has no predecessors or if `n` is 0, the set will contain only the node.

    Example
    _______

    .. code-block:: python

        >>> import networkx as nx
        >>> G = nx.DiGraph()
        >>> G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        >>> get_n_level_ancestors(G, 'D', 2)
        >>> # Output: {'B', 'C', 'D'}

    :notes:
        - The function includes the node itself in the result.
        - Predecessors (ancestors) are determined by traversing the graph from the node upward.
        - If `n` is greater than the number of levels in the graph, the function will return all possible ancestors
          of the node.
    """
    ancestors = set()
    current_level = {node}

    for _ in range(n):
        current_level = {pred for ancestor in current_level for pred in graph.predecessors(ancestor)}
        ancestors.update(current_level)

    ancestors.add(node)
    return ancestors


def get_n_level_descendants(graph: nx.DiGraph, node: str, n: int) -> t.Set[str]:
    """
    Return the set of all descendants of a given node in a directed graph up to `n` levels deep.

    This function traverses `n` levels of the graph (from parent to child), starting from the specified `node`.
    It collects all nodes that are reachable from `node` within `n` levels of successors (descendants).
    The original node is also included in the resulting set.

    :param graph: A directed graph represented by a NetworkX DiGraph object, where nodes represent entities and
        directed edges represent relationships between them.

    :param node: The node for which to find descendants. This is the starting point of the traversal.

    :param n: The number of levels to traverse downwards in the graph, i.e., how far to go in the graph to find
        descendants. If n = 1, it returns immediate successors; if n = 2, it includes successors of successors,
        and so on.

    :return:
        A set containing the node itself and all descendants reachable within `n` levels from the node.
        If the node has no successors or if `n` is 0, the set will contain only the node.

    Example
    _______

    .. code-block:: python

        >>> import networkx as nx
        >>> G = nx.DiGraph()
        >>> G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        >>> get_n_level_descendants(G, 'A', 2)
        >>> # Output: {'A', 'B', 'C'}

    :notes:
        - The function includes the node itself in the result.
        - Successors (descendants) are determined by traversing the graph from the node downward.
        - If `n` is greater than the number of levels in the graph, the function will return all possible descendants
          of the node.
    """
    descendants = set()
    current_level = {node}

    for _ in range(n):
        current_level = {succ for descendant in current_level for succ in graph.successors(descendant)}
        descendants.update(current_level)

    descendants.add(node)
    return descendants


def get_all_ancestors(graph: nx.DiGraph, node: str) -> t.Set[str]:
    """
    Return the set of all ancestors of a given node in a directed graph.

    This function finds all ancestors of the specified `node` by traversing the directed graph upwards.
    It collects all nodes that have a path leading to the `node`. The original node is also included in the resulting set.

    :param graph: A directed graph represented by a NetworkX DiGraph object, where nodes represent entities and
        directed edges represent relationships between them.

    :param node: The node for which to find ancestors. This is the starting point of the traversal.

    :return: A set containing the node itself and all ancestors reachable from the node. If the node has no ancestors,
        the set will contain only the node.

    Example
    _______

    .. code-block:: python

        >>> import networkx as nx
        >>> G = nx.DiGraph()
        >>> G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        >>> get_all_ancestors(G, 'D')
        >>> # Output: {'B', 'C', 'D'}

    :notes:
        - The function includes the node itself in the result.
        - Ancestors are determined by traversing the graph from the node upward.
        - The function returns all ancestors, regardless of how many levels deep they are.
    """
    ancestors = {node}

    for anc in nx.ancestors(graph, node):
        ancestors.add(anc)

    return ancestors


def get_all_descendants(graph: nx.DiGraph, node: str) -> t.Set[str]:
    """
    Return the set of all descendants of a given node in a directed graph.

    This function finds all descendants of the specified `node` by traversing the directed graph downwards.
    It collects all nodes that have a path originating from the `node`. The original node is also included in the
    resulting set.

    :param graph: A directed graph represented by a NetworkX DiGraph object, where nodes represent entities and
        directed edges represent relationships between them.

    :param node: The node for which to find descendants. This is the starting point of the traversal.

    :return: A set containing the node itself and all descendants reachable from the node. If the node has no
        descendants, the set will contain only the node.

    Example
    _______

    .. code-block:: python

        >>> import networkx as nx
        >>> G = nx.DiGraph()
        >>> G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        >>> get_all_descendants(G, 'A')
        >>> # Output: {'A', 'B', 'C', 'D'}

    :notes:
        - The function includes the node itself in the result.
        - Descendants are determined by traversing the graph from the node downward.
        - The function returns all descendants, regardless of how many levels deep they are.
    """
    descendants = {node}

    for desc in nx.descendants(graph, node):
        descendants.add(desc)

    return descendants
