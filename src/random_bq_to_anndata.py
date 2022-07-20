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

# assumes that cas_cell_info.cas_cell_index values are a continuous list of ints
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

# Retrieve a list of all feature objects, ordered by cas_feature_index
def get_features(project, dataset, client):
    sql = f"SELECT cas_feature_index, original_feature_id, feature_name FROM `{project}.{dataset}.cas_feature_info` ORDER BY cas_feature_index"
    query = client.query(sql)
    features = []
    for row in query:
        feature = Feature(row["cas_feature_index"], row["original_feature_id"], row["feature_name"])
        if (len(features) > 0) and (feature.cas_feature_index != features[-1].cas_feature_index + 1):
            # It's not the first element in the list, and it's cas_feature_index is not one more than the previous element's
            raise Exception("ERROR: Non-continuous values for `cas_feature_index` in table `cas_feature_info`")
        features.append(feature)
    return features

# Return a list of cell objects (for the random_cell_ids), ordered by cas_cell_index
def get_cells(project, dataset, client, cell_ids):
    in_clause = f" cas_cell_index IN ({','.join(map(str, cell_ids))})"
    sql = f"SELECT cas_cell_index, original_cell_id, cell_type FROM `{project}.{dataset}.cas_cell_info` WHERE " + in_clause + " ORDER BY cas_cell_index"
    query = client.query(sql)
    cells = []
    for row in query:
        cells.append(Cell(row["cas_cell_index"], row["original_cell_id"], row["cell_type"]))
    return cells

def get_matrix_data(project, dataset, client, random_cell_ids):
    in_clause = f" matrix.cas_cell_index IN ({','.join(map(str, random_cell_ids))})"

    # at some point, we will probably want create temp table of cell_ids and then JOIN on it
    # instead of an IN clause
    # NOTE - really don't need to join to cas_cell_index and cas_feature_index here anymore.
    sql = f"SELECT matrix.cas_cell_index, matrix.cas_feature_index, raw_counts AS count FROM `{project}.{dataset}.cas_cell_info` AS cell, `{project}.{dataset}.cas_raw_count_matrix` AS matrix WHERE matrix.cas_cell_index = cell.cas_cell_index AND" + in_clause + " ORDER BY matrix.cas_cell_index, matrix.cas_feature_index"
    query = client.query(sql)
    return query.result()


def random_bq_to_anndata(project, dataset, num_cells, output_file_prefix):
    client = bigquery.Client(project=project)

    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    random_cell_ids = get_random_cell_ids(project, dataset, client, num_cells)

    # Read the cell information and store for later
    cells = get_cells(project, dataset, client, random_cell_ids)
    original_cell_ids = []
    cell_types = []
    for cell in cells:
        original_cell_ids.append(cell.original_cell_id)
        cell_types.append(cell.cell_type)

    # Read the feature information and store for later.
    features = get_features(project, dataset, client)

    feature_ids = []
    feature_names = []
    for feature in features:
        feature_ids.append(feature.original_feature_id)
        feature_names.append(feature.feature_name)

    # Note that this method requires that the result set returned by get_cell_data be sorted by cas_cell_index (or, could also be sorted by original_cell_id)
    cell_data = get_matrix_data(project, dataset, client, random_cell_ids)

    (index_ptr, indices, data) = generate_sparse_matrix(cell_data, features[0].cas_feature_index, features[-1].cas_feature_index)

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

def generate_sparse_matrix(cell_data, minimum_feature_index, maximum_feature_index):
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

    max_col_num = maximum_feature_index - minimum_feature_index
    max_col_num_populated = False

    last_cas_cell_index = None
    for row in cell_data:
        cas_cell_index = row["cas_cell_index"]
        col_num = row["cas_feature_index"] - minimum_feature_index
        if col_num == max_col_num:
            max_col_num_populated = True

        if cas_cell_index != last_cas_cell_index:
            if last_cas_cell_index is not None:
                # We have just started reading data for a new 'cas_cell_index'.
                index_ptr.append(len(indices))
                print(f"finishing record for cas_cell_index: {last_cas_cell_index}")

            last_cas_cell_index = cas_cell_index

        indices.append(col_num)
        data.append(row["count"])

    # NOTE: In order to satisfy the dimensionality of anndata object that will be generated for this sparse matrix,
    # we need to ensure that there is a value in the sparse data matrix for the MAXIMUM feature index (that is,
    # the right-most possible column. It is possible that there will be NO data in the passed cell_data for this feature,
    # so we check to see if there is any data for that MAXIMUM feature index and if there is not, add a 0.
    if not max_col_num_populated:
        print(f"No data populated for maximum cas_feature_index {maximum_feature_index} - inserting a count of 0 here")
        indices.append(max_col_num)
        data.append(0)

    index_ptr.append(len(indices))
    print(f"finishing record for cas_cell_index: {last_cas_cell_index}")

    return index_ptr, indices, data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)
    parser.add_argument('--output_file_prefix', type=str, help='The prefix of the anndata (.h5ad) file that will be created', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
