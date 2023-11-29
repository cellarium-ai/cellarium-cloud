import typing as t

from google.cloud import bigquery


def _parse_match_object(row: bigquery.Row) -> t.Dict[str, t.Any]:
    """
    Parse a single row of a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param row: Row of the BigQuery job object representing the query execution.

    :return: Dictionary representing the query results.
    """
    return {
        "cell_type": row["cell_type"],
        "cell_count": row["cell_count"],
        "min_distance": row["min_distance"],
        "p25_distance": row["p25_distance"],
        "median_distance": row["median_distance"],
        "p75_distance": row["p75_distance"],
        "max_distance": row["max_distance"],
    }


def _parse_match_object_dev_details(row: bigquery.Row) -> t.Dict[str, t.Any]:
    """
    Parse a single row of a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param row: Row of the BigQuery job object representing the query execution.

    :return: Dictionary representing the query results.
    """
    dataset_ids_with_counts = []
    for dataset_ids_with_counts_struct in row["dataset_ids_with_counts"]:
        dataset_ids_with_counts.append(
            {
                "dataset_id": dataset_ids_with_counts_struct["dataset_id"],
                "count_per_dataset": dataset_ids_with_counts_struct["count_per_dataset"],
                "min_distance": dataset_ids_with_counts_struct["min_distance"],
                "max_distance": dataset_ids_with_counts_struct["max_distance"],
                "median_distance": dataset_ids_with_counts_struct["median_distance"],
                "mean_distance": dataset_ids_with_counts_struct["mean_distance"],
            }
        )
    return {
        "cell_type": row["cell_type"],
        "cell_count": row["cell_count"],
        "min_distance": row["min_distance"],
        "p25_distance": row["p25_distance"],
        "median_distance": row["median_distance"],
        "p75_distance": row["p75_distance"],
        "max_distance": row["max_distance"],
        "dataset_ids_with_counts": dataset_ids_with_counts,
    }


def parse_match_query_job(
    query_job: bigquery.QueryJob, include_dev_details: bool = False
) -> t.List[t.Dict[str, t.Any]]:
    """
    Parse a BigQuery job object representing a query to retrieve metadata for a matching query.

    :param query_job: BigQuery job object representing the query execution.
    :param include_dev_details: Boolean indicating whether to include a breakdown of the number of cells by dataset

    :return: List of dictionaries representing the query results.
    """
    results = []

    last_query_id = None
    data = {}

    for row in query_job:
        if last_query_id is None or last_query_id != row["query_id"]:
            # emit data and reset state if this isn't the first time through
            if last_query_id is not None:
                results.append(data)

            data = {"query_cell_id": row["query_id"], "matches": []}
            last_query_id = data["query_cell_id"]

        if not include_dev_details:
            x = _parse_match_object(row=row)
        else:
            x = _parse_match_object_dev_details(row=row)

        data["matches"].append(x)

    results.append(data)
    return results
