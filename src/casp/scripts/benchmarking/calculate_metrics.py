import concurrent.futures as concurrency
import multiprocessing
import pickle
import traceback
import typing as t
import logging

import anndata
import pandas as pd
import wandb
from smart_open import open

from casp.scripts.benchmarking import utils


logger = logging.getLogger(__name__)


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
    hops = [co_resource[ground_truth_cl_name][f"hop_{i}"] for i in range(num_hops + 1)]
    true_positives = [0.0] * len(hops)
    false_positives = [0.0] * len(hops)

    for match in query_cell_obj["matches"]:
        match_cl_name = match["cell_type_ontology_term_id"]
        match_score = match["score"]
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
    ground_truth_cl_name: str,
    co_resource: t.Dict[str, t.Any],
    num_hops=4,
):
    """
    Calculate performance metrics for a query cell against the ground truth.

    :param query_cell_obj: The query cell object containing ontology aware scores for.
    :param ground_truth_cl_name: The ground truth cell type label.
    :param co_resource: Cell ontology precalculated resource used in CZI.
    :param num_hops: Number of hops to consider.

    :return: A dictionary containing sensitivity, specificity, and F1 score for each hop level.
    """
    if ground_truth_cl_name not in co_resource:
        metrics_na = {
            "query_cell_id": query_cell_obj["query_cell_id"],
            "detail": f"Couldn't find cell type {ground_truth_cl_name} in Cell Ontology resource",
        }
        for hop in range(num_hops + 1):
            for metric in ["sensitivity", "specificity", "f1_score"]:
                metrics_na[f"hop_{hop}_{metric}"] = None

        return metrics_na

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
    ground_truth_cl_names: t.Iterable[str],
    cas_result: t.List[t.Dict[str, t.Any]],
    co_resource: t.Dict[str, t.Any],
    num_hops: int,
) -> pd.DataFrame:
    """
    Calculate performance metrics for CAS (Cell Annotation Service) output.

    This function processes the CAS output and the ground truth cell types from the input AnnData object,
    computes performance metrics (sensitivity, specificity, and F1 score) for each hop level, and returns
    the results as a pandas DataFrame.

    :param ground_truth_cl_names: Iterable containing ground truths to benchmark against
    :param cas_result: The result from the CAS annotation, containing query results and matches for each cell.
    :param co_resource: Cell ontology precalculated resource used in CZI.
    :param num_hops: Number of hops to consider.

    :return: A pandas DataFrame with query cell IDs as the index and sensitivity, specificity, and F1 score
         for each hop level as columns.
    """
    df_result = pd.DataFrame()

    for query_res_obj, ground_truth in zip(cas_result, ground_truth_cl_names):
        query_cell_metrics = calculate_metrics_for_query_cell(
            query_cell_obj=query_res_obj,
            ground_truth_cl_name=ground_truth,
            co_resource=co_resource,
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
    ground_truth_cl_names: t.List[str],
    cas_result: t.Dict[str, t.Any],
    co_resource: t.Dict[str, t.Any],
    num_hops: int,
    batch_size: int = 20000,
) -> pd.DataFrame:
    """
    Calculate metrics for CAS output in batches using multiprocessing.

    :param ground_truth_cl_names: List of ground truth labels.
    :param cas_result: List of CAS result dictionaries.
    :param co_resource: Cell ontology precalculated resource used in CZI.
    :param num_hops: Number of hops for evaluation.
    :param batch_size: Size of each batch. Default is 5000.

    :raises ValueError: If the length of ground_truths and cas_result do not match.

    :return: DataFrame containing the calculated metrics.
    """
    result_dfs = []
    num_workers = multiprocessing.cpu_count()
    with concurrency.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []

        cas_result_batches = split_into_batches(data_list=cas_result["data"], batch_size=batch_size)
        ground_truth_cl_names_batches = split_into_batches(data_list=ground_truth_cl_names, batch_size=batch_size)

        for i, (cas_result_batch, ground_truth_cl_names_batch) in enumerate(
            zip(cas_result_batches, ground_truth_cl_names_batches)
        ):
            calculate_metrics_kwargs = {
                "cas_result": cas_result_batch,
                "ground_truth_cl_names": ground_truth_cl_names_batch,
                "co_resource": co_resource,
                "num_hops": num_hops,
            }
            logger.info(f"Submitting batch {i}. Length: {len(cas_result_batch)}...")
            future = executor.submit(calculate_metrics_for_cas_output, **calculate_metrics_kwargs)
            futures.append(future)

        logger.info("Waiting for metrics to be calculated...")
        done, not_done = concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)
        logger.info("Aggregating results")
        for future in done:
            try:
                # Attempt to get the result of the future
                batch_metrics_df = future.result()
            except Exception as e:
                # If an exception is raised, print the exception details
                logging.error(f"Exception occurred: {e}")
                # Format and print the full traceback
                traceback.print_exception(type(e), e, e.__traceback__)
            else:
                result_dfs.append(batch_metrics_df)

    df_result = pd.concat(result_dfs)
    logger.info("Done.")
    return df_result


def calculate_metrics_for_cas_responses(
    dataset_paths: t.Union[str, t.List[str]],
    cas_result_paths: t.Union[str, t.List[str]],
    model_name: str,
    num_hops: int,
    co_resource_path: str,
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
    :param co_resource_path: Path to the co-resource file for Cell ontology used in CZI.
    :param output_path: Path to save the output metrics.
    :param wandb_project: Weights and Biases project name.
    :param batch_size: Batch size for processing cas results in parallel processes.
    """
    with open(co_resource_path, "rb") as f:
        co_resource = pickle.load(f)

    # Get all columns list for the final result
    metrics = ["sensitivity", "specificity", "f1_score"]
    metric_column_list = [f"hop_{hop}_{metric}" for hop in range(num_hops + 1) for metric in metrics]

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

        # Cell Ontology resource has a different delimiter than what they use in the data.
        ground_truth_cl_names = adata.obs.cell_type_ontology_term_id.map(lambda x: x.replace(":", "_"))

        if len(cas_result.data) != len(ground_truth_cl_names):
            raise ValueError(
                "Length of ground truths array in input dataset does not correspond to the length of cas result list. "
                "Make sure you're using the right dataset and cas_result combination"
            )
        cas_result_dict = cas_result.model_dump()
        df_metrics = calculate_metrics_for_cas_output_in_batches(
            ground_truth_cl_names=ground_truth_cl_names,
            cas_result=cas_result_dict,
            co_resource=co_resource,
            num_hops=num_hops,
            batch_size=batch_size,
        )

        df_metrics.to_csv(metrics_output_filepath)
        results.append(df_metrics)

    means_list = [df[metric_column_list].mean() for df in results]
    means_df = pd.DataFrame(means_list)
    aggregated_means = means_df.mean()
    run.log(aggregated_means.to_dict())
