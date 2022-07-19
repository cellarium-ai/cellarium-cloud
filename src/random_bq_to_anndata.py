from google.cloud import bigquery
import argparse
import random
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

class Cell:
    cas_cell_index = 0
    original_cell_id = ""
    cell_type = ""

    def __init__(self, cas_cell_index, original_cell_id, cell_type):
        self.cas_cell_index = cas_cell_index
        self.original_cell_id = original_cell_id
        self.cell_type = cell_type

class Feature:
    cas_feature_index = 0
    original_feature_id = ""
    feature_name = ""

    def __init__(self, cas_feature_index, original_feature_id, feature_name):
        self.cas_feature_index = cas_feature_index
        self.original_feature_id = original_feature_id
        self.feature_name = feature_name

# assumes that cas_cell_info.cas_cell_index values are a contiguous list of ints
def get_random_cell_ids(project, dataset, client, num_cells):
    query = client.query(f"SELECT MIN(cas_cell_index) AS min_cas_cell_index, MAX(cas_cell_index) AS max_cas_cell_index FROM `{project}.{dataset}.cas_cell_info`")
    row = list(query.result())[0]
    min_cas_cell_index, max_cas_cell_index = row.min_cas_cell_index, row.max_cas_cell_index
    print(f"Getting {num_cells} random IDs between {min_cas_cell_index} and {max_cas_cell_index}...")
    cell_ids = list(range(min_cas_cell_index, max_cas_cell_index + 1))
    random.shuffle(cell_ids)
    del cell_ids[num_cells:]
    print(f"Random IDs: {cell_ids}")
    return cell_ids

# Retrieve all of the features from the database
def get_features(project, dataset, client):
    sql = f"SELECT cas_feature_index, original_feature_id, feature_name FROM `{project}.{dataset}.cas_feature_info` ORDER BY cas_feature_index"
    query = client.query(sql)
    ordered_feature_indices = []
    feature_index_to_feature = {}
    for row in query:
        index = row["cas_feature_index"]
        ordered_feature_indices.append(index)
        feature_index_to_feature[index] = Feature(row["cas_feature_index"], row["original_feature_id"], row["feature_name"])
    return ordered_feature_indices, feature_index_to_feature

# Retrieve metadata for the selected (random) cells from the database
def get_cells(project, dataset, client, random_cell_ids):
    in_clause = f" cas_cell_index IN ({','.join(map(str, random_cell_ids))})"
    sql = f"SELECT cas_cell_index, original_cell_id, cell_type FROM `{project}.{dataset}.cas_cell_info` WHERE " + in_clause + " ORDER BY cas_cell_index"
    query = client.query(sql)
    ordered_cell_indices = []
    cell_index_to_cell = {}
    for row in query:
        index = row["cas_cell_index"]
        ordered_cell_indices.append(index)
        cell_index_to_cell[index] = Cell(row["cas_cell_index"], row["original_cell_id"], row["cell_type"])
    return ordered_cell_indices, cell_index_to_cell

def get_matrix_data(project, dataset, client, random_cell_ids):
    in_clause = f" matrix.cas_cell_index IN ({','.join(map(str, random_cell_ids))})"

    # at some point, we will probably want create temp table of cell_ids and then JOIN on it
    # instead of an IN clause
    # NOTE - really don't need to join to cas_cell_index and cas_feature_index here anymore.
    sql = f"SELECT matrix.cas_cell_index, original_cell_id, matrix.cas_feature_index, original_feature_id, raw_counts AS count FROM `{project}.{dataset}.cas_cell_info` AS cell, `{project}.{dataset}.cas_feature_info` AS feature, `{project}.{dataset}.cas_raw_count_matrix` AS matrix WHERE matrix.cas_cell_index = cell.cas_cell_index AND matrix.cas_feature_index = feature.cas_feature_index AND" + in_clause + " ORDER BY matrix.cas_cell_index, matrix.cas_feature_index"
    query = client.query(sql)
    return query.result()


def random_bq_to_anndata(project, dataset, num_cells, output_file_prefix):
    client = bigquery.Client(project=project)

    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    random_cell_ids = get_random_cell_ids(project, dataset, client, num_cells)

    # Read the cell information and store metadata for later
    ordered_cell_indices, cell_index_to_cell = get_cells(project, dataset, client, random_cell_ids)

    original_cell_ids = []
    cell_types = []
    for cell_index in ordered_cell_indices:
        original_cell_ids.append(cell_index_to_cell[cell_index].original_cell_id)
        cell_types.append(cell_index_to_cell[cell_index].cell_type)

    # Read the feature information and store metadata for later.
    ordered_feature_indices, feature_index_to_feature = get_features(project, dataset, client)

    feature_ids = []
    feature_names = []
    last_feature_index = None
    for feature_index in ordered_feature_indices:
        if (last_feature_index is None):
            last_feature_index = feature_index - 1
        # # We expect all entries in cas_feature_info to be stored in continuous fashion. (1, 2, 3, ...)
        if (feature_index != last_feature_index + 1):
            raise Exception("ERROR: Non-continuous values for `cas_feature_index` in table `cas_feature_info`")

        feature_ids.append(feature_index_to_feature[feature_index].original_feature_id)
        feature_names.append(feature_index_to_feature[feature_index].feature_name)
        last_feature_index = feature_index

    # Note that this method requires that the result set returned by get_cell_data be sorted by cas_cell_index (or, could also be sorted by original_cell_id)
    cell_data = get_matrix_data(project, dataset, client, random_cell_ids)

    (index_ptr, indices, data) = generate_sparse_matrix(cell_data, ordered_feature_indices[0])

    # Create the matrix from the sparse data representation generated above.
    counts = csr_matrix((data, indices, index_ptr), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = original_cell_ids
    adata.obs["cell_type"] = cell_types
    adata.var_names = feature_ids
    adata.var["feature_name"] = feature_names

    # See https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html?highlight=raw#anndata.AnnData.raw
    #  for why we set 'raw' thusly: "The raw attribute is initialized with the current content of an object by setting:"
    adata.raw = adata
    adata.write(f'{output_file_prefix}.h5ad', compression="gzip")

def generate_sparse_matrix(cell_data, minimum_feature_index):
    # For representation of the sparse matrix
    index_ptr = [0]
    indices = []
    data = []

    # Builds the sparse matrix for storing the original_cell_ids (as obs) and features (as var).
    # From: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
    # If we want to store the following sparse matrix:
    # OBS \ VAR
    #       ENS1 ENS2 ENS3 ENS4
    # AAA   1    7     2   n/a
    # AAC   0    n/a   4   n/a
    # AAG   3    n/a   0   1
    #
    # we represent it as:
    # index_ptr = [0, 3, 5, 8]
    # indices   = [0, 1, 2, 0, 2, 0, 2, 3]
    # data      = [1, 7, 2, 0, 4, 3, 0, 1]

    last_original_cell_id = None
    for row in cell_data:
        original_cell_id = row["original_cell_id"]
        col_num = row["cas_feature_index"] - minimum_feature_index
        count = row["count"]

        if original_cell_id != last_original_cell_id:
            if last_original_cell_id is not None:
                # We have just started reading data for a new 'original_cell_id'.
                index_ptr.append(len(indices))
                # print(f"finishing record for original_cell_id: {last_original_cell_id}")
                # print(f"index_ptr: {index_ptr}")
            last_original_cell_id = original_cell_id

        indices.append(col_num)
        data.append(count)

    index_ptr.append(len(indices))
    # print(f"finishing record for original_cell_id: {last_original_cell_id}")
    # print(f"index_ptr: {index_ptr}")

    return index_ptr, indices, data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)
    parser.add_argument('--output_file_prefix', type=str, help='The prefix of the anndata (.h5ad) file that will be created', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
