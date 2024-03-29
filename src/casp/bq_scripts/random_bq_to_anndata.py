"""
Selects a random subset of cells of a specified size from BigQuery and writes the data to output AnnData files.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
from google.cloud import bigquery
from scipy.sparse import coo_matrix


class Cell:
    """
    Represents a CASP cell.
    """

    def __init__(self, cas_cell_index, original_cell_id, cell_type, obs_metadata, cas_ingest_id):
        self.cas_cell_index = cas_cell_index
        self.original_cell_id = original_cell_id
        self.cell_type = cell_type
        self.obs_metadata = obs_metadata
        self.cas_ingest_id = cas_ingest_id


class Feature:
    """
    Represents a CASP feature / gene.
    """

    def __init__(self, cas_feature_index, original_feature_id, feature_name, var_metadata):
        self.cas_feature_index = cas_feature_index
        self.original_feature_id = original_feature_id
        self.feature_name = feature_name
        self.var_metadata = var_metadata


# Returns an array of <num_cells> cells pulled at random from the cas_cell_info table, sorted by `cas_cell_index`.
def get_random_cells(project, dataset, client, num_cells):
    """
    Selects the specified number of cells randomly from the specified dataset.
    """
    query = f"""

        SELECT cas_cell_index, original_cell_id, cell_type, obs_metadata, cas_ingest_id, rand() AS rand_val FROM
            `{project}.{dataset}.cas_cell_info` ORDER BY rand_val LIMIT {num_cells}

    """

    result = client.query(query)
    cells = [
        Cell(
            cas_cell_index=row["cas_cell_index"],
            original_cell_id=row["original_cell_id"],
            cell_type=row["cell_type"],
            obs_metadata=row["obs_metadata"],
            cas_ingest_id=row["cas_ingest_id"],
        )
        for row in result
    ]
    cells.sort(key=lambda c: c.cas_cell_index)
    return cells


# Retrieve a list of all feature objects, ordered by cas_feature_index
def get_features(project, dataset, client, ingest_id):
    """
    Retrieve all features in the dataset for the specified ingest_id.
    """
    sql = f"""

    SELECT cas_feature_index, original_feature_id, feature_name, var_metadata FROM
        `{project}.{dataset}.cas_feature_info` WHERE cas_ingest_id = "{ingest_id}" ORDER BY cas_feature_index

    """

    query = client.query(sql)
    features = []
    for row in query:
        features.append(
            Feature(row["cas_feature_index"], row["original_feature_id"], row["feature_name"], row["var_metadata"])
        )
    return features


def get_ingest_info(project, dataset, client, ingest_id):
    """
    Retrieve ingest info for the specified ingest ids.
    """
    sql = f"""

    SELECT cas_ingest_id, uns_metadata FROM `{project}.{dataset}.cas_ingest_info` WHERE
        cas_ingest_id = "{ingest_id}"

    """

    for row in client.query(sql):
        return row["uns_metadata"]


def get_matrix_data(project, dataset, client, cells):
    """
    Retrieve raw matrix data for the specified cell ids.
    """
    str_cell_ids = map(str, [c.cas_cell_index for c in cells])
    in_clause = f" cas_cell_index IN ({','.join(str_cell_ids)})"

    # An IN clause might seem like it wouldn't scale well but in practice seems to perform okay to at least ~100 K
    # random cells. Overall 100 K performance is poor but this doesn't seem to be a major contributor to that.
    sql = f"""

    SELECT cas_cell_index, cas_feature_index, raw_counts AS count FROM
        `{project}.{dataset}.cas_raw_count_matrix` WHERE {in_clause} ORDER BY cas_cell_index, cas_feature_index

    """
    query = client.query(sql)
    return query.result()


def assign_obs_var_metadata(dataframe, json_strings):
    """
    :param dataframe: Pandas DataFrame into which metadata should be written back
    :param json_strings: An iterable of stringified JSON to be written back to the DataFrame
    """
    jsons = [json.loads(s) for s in json_strings]

    # The keys will be the same for all JSONs so iterate over the keys of the first:
    for key in jsons[0].keys():
        values = [j[key] for j in jsons]
        dataframe[key] = values


def assign_uns_metadata(anndata, json_string):
    """
    Assign the specified JSON string into the `uns` metadata for this AnnData file.
    """
    json_dict = json.loads(json_string)
    for key, val in json_dict.items():
        anndata.uns[key] = val


def random_bq_to_anndata(project, dataset, num_cells, output_file_prefix):
    """
    Main function to select a specified number of random cells from the BigQuery dataset and write the associated data
    to an AnnData file.
    """
    client = bigquery.Client(project=project)

    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    random_cells = get_random_cells(project, dataset, client, num_cells)

    cells_by_ingest_id = defaultdict(list)
    for cell in random_cells:
        cells_by_ingest_id[cell.cas_ingest_id].append(cell)

    for ingest_id, cells in cells_by_ingest_id.items():
        id_only = ingest_id.split("-")[-1]
        output_file_name = f"{output_file_prefix}-{id_only}.h5ad"
        print(f"Writing {len(cells)} cells' data to '{output_file_name}'...")

        original_cell_ids = []
        cas_cell_index_to_row_num = {}
        for row_num, cell in enumerate(cells):
            original_cell_ids.append(cell.original_cell_id)
            cas_cell_index_to_row_num[cell.cas_cell_index] = row_num

        # Read the feature information and store for later.
        features = get_features(project, dataset, client, ingest_id)

        feature_ids = []
        cas_feature_index_to_col_num = {}
        for col_num, feature in enumerate(features):
            feature_ids.append(feature.original_feature_id)
            cas_feature_index_to_col_num[feature.cas_feature_index] = col_num

        matrix_data = get_matrix_data(project, dataset, client, cells)

        rows, columns, data = convert_matrix_data_to_coo_matrix_input_format(
            matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num
        )
        # Create the matrix from the sparse data representation generated above.
        counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(features)), dtype=np.float32)

        # Convert the COO matrix to CSR for loading into AnnData
        adata = ad.AnnData(counts.tocsr(copy=False))
        adata.obs.index = original_cell_ids
        adata.var.index = feature_ids
        assign_obs_var_metadata(adata.obs, [c.obs_metadata for c in cells])
        assign_obs_var_metadata(adata.var, [f.var_metadata for f in features])

        ingest_metadata = get_ingest_info(project, dataset, client, ingest_id)
        assign_uns_metadata(adata, ingest_metadata)

        # See https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html?highlight=raw#anndata.AnnData.raw
        # for why we set 'raw' thusly: "The raw attribute is initialized with the current content of an object by setting:"
        adata.raw = adata
        adata.write(Path(output_file_name), compression="gzip")
    print("Done.")


def convert_matrix_data_to_coo_matrix_input_format(
    matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num
):
    """
    Convert raw matrix data to coordinate format suitable for writing out an AnnData file.
    """
    rows = []
    columns = []
    data = []

    for row in matrix_data:
        cas_cell_index = row["cas_cell_index"]
        cas_feature_index = row["cas_feature_index"]
        try:
            row_num = cas_cell_index_to_row_num[cas_cell_index]
        except KeyError as exc:
            raise Exception(
                f"ERROR: Unable to find entry for cas_cell_index: {cas_cell_index} in lookup table"
            ) from exc

        try:
            col_num = cas_feature_index_to_col_num[cas_feature_index]
        except KeyError as exc:
            raise Exception(
                f"ERROR: Unable to find entry for cas_feature_index: {cas_feature_index} in lookup table"
            ) from exc

        rows.append(row_num)
        columns.append(col_num)
        data.append(row["count"])

    return rows, columns, data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, description="Query CASP tables for random cells")
    parser.add_argument("--project", type=str, help="BigQuery Project", required=True)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--num_cells", type=int, help="Number of cells to return", required=True)
    parser.add_argument(
        "--output_file_prefix",
        type=str,
        help="The prefix of the anndata (.h5ad) file that will be created",
        required=True,
    )

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
