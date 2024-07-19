import concurrent.futures as concurrency
import multiprocessing
import pickle
import traceback
import typing as t

import anndata
import pandas as pd
import wandb
from smart_open import open

from casp.scripts.benchmarking import utils


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
    if ground_truth not in co_resource:
        metrics_na = {
            "query_cell_id": query_cell_obj["query_cell_id"],
            "detail": f"Couldn't find cell type {ground_truth} in Cell Ontology resource"
        }
        for i in range(num_hops):
            metrics_na[f"hop_{i}_sensitivity"] = None
            metrics_na[f"hop_{i}_specificity"] = None
            metrics_na[f"hop_{i}_f1_score"] = None

        return metrics_na

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

    query_cell_metrics = {"query_cell_id": query_cell_obj["query_cell_id"], "detail": ""}

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
        query_cell_metrics = calculate_metrics_for_query_cell(
            query_cell_obj=query_res_obj,
            ground_truth=ground_truth,
            co_resource=co_resource,
            cl_labels_to_names_map=cl_labels_to_names_map,
            num_hops=num_hops,
        )

        df_result = pd.concat([df_result, pd.DataFrame([query_cell_metrics])])

    df_result = df_result.set_index("query_cell_id")
    return df_result


def split_into_batches(data_list, batch_size):
    """
    Split a list into batches of a given size.

    :param data_list: List of data to be split into batches.
    :param batch_size: The size of each batch.

    :return: List of batches, where each batch is a list.
    """
    return [data_list[i : i + batch_size] for i in range(0, len(data_list), batch_size)]


def calculate_metrics_for_cas_output_in_batches(
    ground_truths: t.List[str],
    cas_result: t.List[t.Dict[str, t.Any]],
    co_resource: t.Dict[str, t.Any],
    cl_labels_to_names_map: t.Dict[str, str],
    num_hops: int,
    batch_size: int = 20000,
) -> pd.DataFrame:
    """
    Calculate metrics for CAS output in batches using multiprocessing.

    :param ground_truths: List of ground truth labels.
    :param cas_result: List of CAS result dictionaries.
    :param co_resource: Cell Ontology precalculated resource dictionary.
    :param cl_labels_to_names_map: Mapping from class labels to cell ontology names.
    :param num_hops: Number of hops for evaluation.
    :param batch_size: Size of each batch. Default is 5000.

    :raises ValueError: If the length of ground_truths and cas_result do not match.

    :return: DataFrame containing the calculated metrics.
    """
    if len(cas_result) != len(ground_truths):
        raise ValueError(
            "Length of ground truths array does not correspond to the length of cas result list. "
            "Make sure you're using the right dataset's ground truths and cas_result combination"
        )
    result_dfs = []
    num_workers = multiprocessing.cpu_count()
    with concurrency.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []

        cas_result_batches = split_into_batches(data_list=cas_result, batch_size=batch_size)
        ground_truths_batches = split_into_batches(data_list=ground_truths, batch_size=batch_size)

        for i, (cas_result_batch, ground_truths_batch) in enumerate(zip(cas_result_batches, ground_truths_batches)):
            calculate_metrics_kwargs = {
                "cas_result": cas_result_batch,
                "ground_truths": ground_truths_batch,
                "co_resource": co_resource,
                "cl_labels_to_names_map": cl_labels_to_names_map,
                "num_hops": num_hops,
            }
            print(f"Submitting batch {i}. Length: {len(cas_result_batch)}...")
            future = executor.submit(calculate_metrics_for_cas_output, **calculate_metrics_kwargs)
            futures.append(future)

        print("Waiting for metrics to be calculated...")
        done, not_done = concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)
        print("Aggregating results")
        for future in done:
            try:
                # Attempt to get the result of the future
                batch_metrics_df = future.result()
            except Exception as e:
                # If an exception is raised, print the exception details
                print(f"Future: {future}")
                print(f"Exception type: {type(e).__name__}")
                print(f"Exception message: {e}")
                # Format and print the full traceback
                traceback.print_exception(type(e), e, e.__traceback__)
            else:
                result_dfs.append(batch_metrics_df)

    df_result = pd.concat(result_dfs)
    print("Done.")
    return df_result


def calculate_metrics_for_cas_responses(
    dataset_paths: t.Union[str, t.List[str]],
    cas_result_paths: t.Union[str, t.List[str]],
    model_name: str,
    num_hops: int,
    co_resource_path: str,
    cl_labels_to_names_map_path: str,
    output_path: str,
    wandb_project: str,
    batch_size: int,
) -> None:
    """
    Calculate metrics for CAS responses and log results to Weights and Biases.

    :param dataset_paths: Paths to the dataset files. Can be a list of paths or a path to txt file with paths.
    :param cas_result_paths: Paths to the CAS result files. Can be a list of paths or a path to txt file with paths.
    :param model_name: Name of the model being evaluated.
    :param num_hops: Number of hops for evaluation.
    :param co_resource_path: Path to the co-resource file.
    :param cl_labels_to_names_map_path: Path to the class labels to names map file.
    :param output_path: Path to save the output metrics.
    :param wandb_project: Weights and Biases project name.
    :param batch_size: Batch size for processing cas results in parallel processes.
    """
    with open(co_resource_path, "rb") as f:
        co_resource = pickle.load(f)

    with open(cl_labels_to_names_map_path, "rb") as f:
        cl_labels_to_names_map = pickle.load(f)

    _dataset_paths = utils.get_paths(paths=dataset_paths)
    _cas_result_paths = utils.get_paths(paths=cas_result_paths)

    run = wandb.init(project=wandb_project, name=f"benchmarking_model_{model_name}")

    metrics_dir = f"{output_path}/{model_name}"
    results = []

    for dataset_file_path, cas_result_path in zip(_dataset_paths, _cas_result_paths):
        dataset_file_name = dataset_file_path.split("/")[-1].split(".")[0]
        metrics_output_filepath = f"{metrics_dir}/metrics_{dataset_file_name}.csv"

        with open(dataset_file_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        with open(cas_result_path, "rb") as f:
            cas_result = pickle.load(f)

        if adata.raw is not None:
            # Getting raw expression counts instead of normalized.
            adata = adata.raw.to_adata()

        ground_truths = adata.obs.cell_type.values
        if len(cas_result) != len(ground_truths):
            raise ValueError(
                "Length of ground truths array in input dataset does not correspond to the length of cas result list. "
                "Make sure you're using the right dataset and cas_result combination"
            )
        df_metrics = calculate_metrics_for_cas_output_in_batches(
            ground_truths=ground_truths,
            cas_result=cas_result,
            co_resource=co_resource,
            cl_labels_to_names_map=cl_labels_to_names_map,
            num_hops=num_hops,
            batch_size=batch_size,
        )

        df_metrics.to_csv(metrics_output_filepath)
        results.append(df_metrics)

    means_list = [df.mean() for df in results]
    means_df = pd.DataFrame(means_list)
    aggregated_means = means_df.mean()
    run.log(aggregated_means.to_dict())
