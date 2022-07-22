import anndata as ad
import argparse
from google.cloud import bigquery
import json
import numpy as np
from pathlib import Path
from scipy.sparse import coo_matrix


class Cell:
    def __init__(self, cas_cell_index, original_cell_id, cell_type, obs_metadata, ingest_id):
        self.cas_cell_index = cas_cell_index
        self.original_cell_id = original_cell_id
        self.cell_type = cell_type
        self.obs_metadata = obs_metadata
        self.ingest_id = ingest_id


class Feature:
    def __init__(self, cas_feature_index, original_feature_id, feature_name, var_metadata):
        self.cas_feature_index = cas_feature_index
        self.original_feature_id = original_feature_id
        self.feature_name = feature_name
        self.var_metadata = var_metadata


# Returns an array of <num_cells> cas_cell_index values pulled at random from the cas_cell_info table
def get_random_cell_ids(project, dataset, client, num_cells):
    query = client.query("SELECT cas_cell_index, rand() AS rand_val FROM " +
                         f"`{project}.{dataset}.cas_cell_info` ORDER BY rand_val LIMIT {num_cells}")
    cas_cell_indices = []
    for row in query:
        cas_cell_indices.append(row["cas_cell_index"])
    print(f"Random IDs: {cas_cell_indices}")
    return cas_cell_indices


# Retrieve a list of all feature objects, ordered by cas_feature_index
def get_features(project, dataset, client):
    sql = "SELECT cas_feature_index, original_feature_id, feature_name, var_metadata FROM " + \
          f"`{project}.{dataset}.cas_feature_info` ORDER BY cas_feature_index"
    query = client.query(sql)
    features = []
    for row in query:
        features.append(Feature(row["cas_feature_index"],
                                row["original_feature_id"],
                                row["feature_name"],
                                row["var_metadata"]))
    return features


# Given a list of cell_ids (= cas_cell_index) this method will return a list of cell objects for them, ordered by
# cas_cell_index
def get_cells(project, dataset, client, cell_ids):
    in_clause = f" cas_cell_index IN ({','.join(map(str, cell_ids))})"
    sql = "SELECT cas_cell_index, original_cell_id, cell_type, obs_metadata, cas_ingest_id FROM " + \
          f"`{project}.{dataset}.cas_cell_info` WHERE {in_clause} ORDER BY cas_cell_index"
    query = client.query(sql)
    cells = []
    for row in query:
        cells.append(Cell(row["cas_cell_index"],
                          row["original_cell_id"],
                          row["cell_type"],
                          row["obs_metadata"],
                          row["cas_ingest_id"]))
    return cells


def get_ingest_info(project, dataset, client, ingest_ids):
    ids = ', '.join([f'"{i}"' for i in ingest_ids])
    sql = f"""SELECT cas_ingest_id, uns_metadata FROM `{project}.{dataset}.cas_ingest_info` WHERE """ + \
          f"""cas_ingest_id IN ({ids})"""

    info = {}
    for row in client.query(sql):
        info[row['cas_ingest_id']] = row['uns_metadata']
    return info


def get_matrix_data(project, dataset, client, random_cell_ids):
    in_clause = f" cas_cell_index IN ({','.join(map(str, random_cell_ids))})"

    # If the number of cells to be selected becomes large, we may want to create a temp table of the cell_ids and then
    # JOIN on it instead of using an IN clause.
    sql = "SELECT cas_cell_index, cas_feature_index, raw_counts AS count FROM " + \
          f"`{project}.{dataset}.cas_raw_count_matrix` WHERE {in_clause} ORDER BY cas_cell_index, cas_feature_index"
    query = client.query(sql)
    return query.result()


def assign_obs_var_metadata(dataframe, json_strings):
    """
    :param dataframe: Pandas DataFrame into which metadata should be written back
    :param json_strings: An iterable of strings stringified JSON to be written back to the DataFrame
    :return:
    """
    jsons = [json.loads(s) for s in json_strings]

    # The keys will be the same for all JSONs so iterate over the keys of the first:
    for key in jsons[0].keys():
        values = [j[key] for j in jsons]
        dataframe[key] = values


def assign_uns_metadata(anndata, json_string):
    json_dict = json.loads(json_string)
    for k, v in json_dict.items():
        anndata.uns[k] = v


def random_bq_to_anndata(project, dataset, num_cells, output_file_prefix):
    client = bigquery.Client(project=project)

    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    random_cell_ids = get_random_cell_ids(project, dataset, client, num_cells)

    # Read the cell information and store for later.
    cells = get_cells(project, dataset, client, random_cell_ids)
    original_cell_ids = []
    cas_cell_index_to_row_num = {}
    for row_num, cell in enumerate(cells):
        original_cell_ids.append(cell.original_cell_id)
        cas_cell_index_to_row_num[cell.cas_cell_index] = row_num

    # Read the feature information and store for later.
    features = get_features(project, dataset, client)

    feature_ids = []
    cas_feature_index_to_col_num = {}
    for col_num, feature in enumerate(features):
        feature_ids.append(feature.original_feature_id)
        cas_feature_index_to_col_num[feature.cas_feature_index] = col_num

    # Note that this method requires that the result set returned by get_cell_data be sorted by cas_cell_index
    matrix_data = get_matrix_data(project, dataset, client, random_cell_ids)

    rows, columns, data = convert_matrix_data_to_coo_matrix_input_format(matrix_data, cas_cell_index_to_row_num,
                                                                         cas_feature_index_to_col_num)
    # Create the matrix from the sparse data representation generated above.
    counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(features)), dtype=np.float32)

    # Convert the COO matrix to CSR for loading into AnnData
    adata = ad.AnnData(counts.tocsr(copy=False))
    adata.obs.index = original_cell_ids
    adata.var.index = feature_ids
    assign_obs_var_metadata(adata.obs, [c.obs_metadata for c in cells])
    assign_obs_var_metadata(adata.var, [f.var_metadata for f in features])
    # TODO make this work with random subsets that span multiple ingested AnnData files
    ingest_metadata = get_ingest_info(project, dataset, client, [cells[0].ingest_id])
    assign_uns_metadata(adata, ingest_metadata[cells[0].ingest_id])

    # See https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html?highlight=raw#anndata.AnnData.raw
    # for why we set 'raw' thusly: "The raw attribute is initialized with the current content of an object by setting:"
    adata.raw = adata
    adata.write(Path(f'{output_file_prefix}.h5ad'), compression="gzip")


def convert_matrix_data_to_coo_matrix_input_format(cell_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num):
    rows = []
    columns = []
    data = []

    for row in cell_data:
        cas_cell_index = row["cas_cell_index"]
        cas_feature_index = row["cas_feature_index"]
        try:
            row_num = cas_cell_index_to_row_num[cas_cell_index]
        except KeyError:
            raise Exception(f"ERROR: Unable to find entry for cas_cell_index: {cas_cell_index} in lookup table")

        try:
            col_num = cas_feature_index_to_col_num[cas_feature_index]
        except KeyError:
            raise Exception(f"ERROR: Unable to find entry for cas_feature_index: {cas_feature_index} in lookup table")

        rows.append(row_num)
        columns.append(col_num)
        data.append(row["count"])

    return rows, columns, data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)
    parser.add_argument('--output_file_prefix', type=str,
                        help='The prefix of the anndata (.h5ad) file that will be created', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
