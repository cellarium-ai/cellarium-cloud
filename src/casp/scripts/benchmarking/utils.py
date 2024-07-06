import typing as t

import pandas as pd


def create_hop_data(node_cl_name, level, co_resource) -> t.Dict[str, t.Any]:
    """
    Create Hop cell ontology data based on `co_resource`. The output is a dictionary with the following keys:

    * `nodes` – A set of all nodes in the current hop.
    * `all_ancestors` – A set of all the ancestors of all the nodes in the current hop.
    * `all_descendants` – A set of all the descendants of all the nodes in the current hop.

    :param node_cl_name: The current node used as the starting point for the hops. This node will always be the only
        node in `nodes` if the level is 0.
    :param level: The level defines the distance used to determine how far we want to traverse in the Cell Ontology
        graph from `node_cl_name`.
    :param co_resource: Cell Ontology helper resource which has already precalculated top_n_level descendants and
        ancestors.
    """
    current_node_co_data = co_resource[node_cl_name]
    nodes = current_node_co_data[f"top_{level}_ancestors"].union(current_node_co_data[f"top_{level}_descendants"])
    all_ancestors = set().union(*[co_resource[hop_cl_name]["all_ancestors"] for hop_cl_name in nodes])
    all_descendants = set().union(*[co_resource[hop_cl_name]["all_descendants"] for hop_cl_name in nodes])

    return {"nodes": nodes, "all_ancestors": all_ancestors, "all_descendants": all_descendants}


def calculate_precision(tp: float, fp: float) -> float:
    """
    Calculate precision.

    :param tp: True positives.
    :param fp: False positives.
    :return: Precision value.
    """
    try:
        return tp / (tp + fp)
    except ZeroDivisionError:
        return 0.0


def calculate_f1(precision: float, recall: float) -> float:
    """
    Calculate F1 score.

    :param precision: Precision value.
    :param recall: Recall value.

    :return: F1 score.
    """
    try:
        return (2 * precision * recall) / (precision + recall)
    except ZeroDivisionError:
        return 0.0


def calculate_tps_and_fps(
    query_cell_obj: t.Dict[str, t.Any], ground_truth_cl_name: str, num_hops: int, co_resource: t.Dict[str, t.Any]
) -> t.Tuple[t.List[float], t.List[float]]:
    """
    Calculate true positives and false positives for each hop level.

    :param query_cell_obj: The query cell object containing ontology aware scores for.
    :param ground_truth_cl_name: The ground truth cell type name.
    :param num_hops: Number of hops to consider.
    :param co_resource: Cell Ontology helper resource which has already precalculated top_n_level descendants and
        ancestors.

    :return: A tuple of lists containing true positive and false positive scores for each hop level.
    """
    hops = [
        create_hop_data(node_cl_name=ground_truth_cl_name, level=i, co_resource=co_resource)
        for i in range(num_hops + 1)
    ]
    true_positives = [0.0] * len(hops)
    false_positives = [0.0] * len(hops)

    for match in query_cell_obj["matches"]:
        match_cl_name = match["cell_type_ontology_term_id"]
        match_score = match["score"]  # Ensure match_score is retrieved here
        match_co_data = co_resource[match_cl_name]
        match_ancestors = match_co_data["all_ancestors"]

        for i, hop in enumerate(hops):
            hop_match_intersect = hop["nodes"].intersection(match_ancestors)
            hop_all_descendants = hop["all_descendants"]
            hop_all_ancestors = hop["all_ancestors"]

            if match_cl_name in hop_match_intersect:
                true_positives[i] = max(match_score, true_positives[i])
            elif match_cl_name not in hop_all_descendants.union(hop_all_ancestors):
                false_positives[i] = max(match_score, false_positives[i])

    return true_positives, false_positives


def calculate_metrics_for_query_cell(
    query_cell_obj: t.Dict[str, t.Any],
    ground_truth: str,
    co_resource: t.Dict[str, t.Set],
    cl_labels_to_names_map: t.Dict[str, str],
    num_hops=4,
):
    """
    Calculate performance metrics for a query cell against the ground truth.

    :param query_cell_obj: The query cell object containing ontology aware scores for.
    :param ground_truth: The ground truth cell type label.
    :param co_resource: Cell Ontology helper resource which has already precalculated top_n_level descendants and
        ancestors.
    :param cl_labels_to_names_map: A dictionary resource that maps cell type names to cl labels
    :param num_hops: Number of hops to consider.

    :return: A dictionary containing sensitivity, specificity, and F1 score for each hop level.
    """
    ground_truth_cl_name = cl_labels_to_names_map[ground_truth]

    true_positives, false_positives = calculate_tps_and_fps(
        query_cell_obj=query_cell_obj,
        ground_truth_cl_name=ground_truth_cl_name,
        num_hops=num_hops,
        co_resource=co_resource,
    )

    sensitivities = [tp for tp in true_positives]
    specificities = [1 - fp for fp in false_positives]
    precisions = [calculate_precision(tp=tp, fp=fp) for tp, fp, in zip(true_positives, false_positives)]
    f1_scores = [
        calculate_f1(precision=precision, recall=sensitivity)
        for precision, sensitivity in zip(precisions, sensitivities)
    ]

    query_cell_metrics = {}

    for i, (sensitivity, specificity, f1_score) in enumerate(zip(sensitivities, specificities, f1_scores)):
        query_cell_metrics[f"hop_{i}_sensitivity"] = sensitivity
        query_cell_metrics[f"hop_{i}_specificity"] = specificity
        query_cell_metrics[f"hop_{i}_f1_score"] = f1_score

    return query_cell_metrics


def calculate_metrics_for_cas_output(
    ground_truths: t.Iterable[str],
    cas_result: t.List[t.Dict[str, t.Any]],
    co_resource: t.Dict[str, t.Any],
    cl_labels_to_names_map: t.Dict[str, str],
    num_hops: int,
) -> pd.DataFrame:
    """
    Calculate performance metrics for CAS (Cell Annotation Service) output.

    This function processes the CAS output and the ground truth cell types from the input AnnData object,
    computes performance metrics (sensitivity, specificity, and F1 score) for each hop level, and returns
    the results as a pandas DataFrame.

    :param ground_truths: Iterable containing ground truths to benchmark against
    :param cas_result: The result from the CAS annotation, containing query results and matches for each cell.
    :param co_resource: Cell Ontology helper resource which has already precalculated top_n_level descendants and
        ancestors.
    :param cl_labels_to_names_map: A dictionary resource that maps cell type names to cl labels
    :param num_hops: Number of hops to consider.

    :return: A pandas DataFrame with query cell IDs as the index and sensitivity, specificity, and F1 score
         for each hop level as columns.
    """
    df_result = pd.DataFrame()

    for query_res_obj, ground_truth in zip(cas_result, ground_truths):

        metrics = calculate_metrics_for_query_cell(
            query_cell_obj=query_res_obj,
            ground_truth=ground_truth,
            co_resource=co_resource,
            cl_labels_to_names_map=cl_labels_to_names_map,
            num_hops=num_hops,
        )

        df_row = {"query_cell_id": query_res_obj["query_cell_id"], **metrics}

        df_result = pd.concat([df_result, pd.DataFrame([df_row])])

    df_result = df_result.set_index("query_cell_id")
    return df_result
