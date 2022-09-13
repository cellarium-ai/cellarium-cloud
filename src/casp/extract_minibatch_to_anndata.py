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

import datetime

from google.cloud.bigquery_storage import BigQueryReadClient
from google.cloud.bigquery_storage import types

class Cell:
    """
    Represents a CASP cell.
    """

    def __init__(self, cas_cell_index):
        self.cas_cell_index = cas_cell_index


class Feature:
    """
    Represents a CASP feature / gene.
    """

    def __init__(self, cas_feature_index, feature_name):
        self.cas_feature_index = cas_feature_index
        self.feature_name = feature_name



# Retrieve a list of all feature objects, ordered by cas_feature_index
def get_features(project, dataset, extract_feature_table, client):
    """
    Retrieve all features in the dataset for this extract
    """
    sql = f"""

    SELECT cas_feature_index, feature_name FROM
        `{project}.{dataset}.{extract_feature_table}` ORDER BY cas_feature_index

    """

    query = client.query(sql)
    features = []
    for row in query:
        features.append(
            Feature(row["cas_feature_index"], row["feature_name"])
        )
    return features

# Retrieve a list of all cells in the bin, ordered by cas_cell_index
def get_cells(project, dataset, extract_cell_table, start_bin, end_bin, client):
    """
    Retrieve all cells in the dataset for this extract bin
    """
    sql = f"""

    SELECT cas_cell_index FROM
        `{project}.{dataset}.{extract_cell_table}` WHERE extract_bin BETWEEN {start_bin} AND {end_bin} ORDER BY cas_cell_index

    """

    result = client.query(sql)
    cells = [
        Cell(
            cas_cell_index=row["cas_cell_index"]
        )
        for row in result
    ]
    return cells
        
def get_matrix_data(project, dataset, table, start_bin, end_bin, client):
    client = BigQueryReadClient()
    table = f"projects/{project}/datasets/{dataset}/tables/{table}"
                
    requested_session = types.ReadSession()
    requested_session.table = table
    # This API can also deliver data serialized in Apache Arrow format.
    # This example leverages Apache Avro.
    requested_session.data_format = types.DataFormat.AVRO
    
    requested_session.read_options.selected_fields = ["cas_cell_index", "feature_data"]
    requested_session.read_options.row_restriction = f"extract_bin BETWEEN {start_bin} AND {end_bin}"
    
    parent = f"projects/{project}"
    session = client.create_read_session(
        parent=parent,
        read_session=requested_session,
        # We'll use only a single stream for reading data from the table. However,
        # if you wanted to fan out multiple readers you could do so by having a
        # reader process each individual stream.
        max_stream_count=1,
    )
    
    print(f"Estimated Bytes to scan {session.estimated_total_bytes_scanned}")
    reader = client.read_rows(session.streams[0].name)
    
    rows = reader.rows(session)
    return rows

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


def extract_minibatch_to_anndata(project, dataset, extract_table_prefix, start_bin, end_bin, output):
    """
    Main function to extract a minibatch from a prepared training extract and write the associated data
    to an AnnData file.
    """
    client = bigquery.Client(project=project)

    print(f"Getting extract bin {start_bin}-{end_bin} data from {project}.{dataset}.{extract_table_prefix}*...")

    # Read the feature information and store for later.
    print(f"Extracting Feature Info...")
    features = get_features(project, dataset, f"{extract_table_prefix}__extract_feature_info", client)

    feature_ids = []
    cas_feature_index_to_col_num = {}
    for col_num, feature in enumerate(features):
        feature_ids.append(feature.feature_name)
        cas_feature_index_to_col_num[feature.cas_feature_index] = col_num


    print(f"Extracting Cell Info...")
    cells = get_cells(project, dataset, f"{extract_table_prefix}__extract_cell_info", start_bin, end_bin, client)

    original_cell_ids = []
    cas_cell_index_to_row_num = {}
    for row_num, cell in enumerate(cells):
        original_cell_ids.append(cell.cas_cell_index)        
        cas_cell_index_to_row_num[cell.cas_cell_index] = row_num
        
    print(f"Extracting Matrix Data...")
    
    matrix_data = get_matrix_data(project, dataset, f"{extract_table_prefix}__extract_raw_count_matrix", start_bin, end_bin, client)

    print(f"Converting Matrix Data to COO format...")

    rows, columns, data = convert_matrix_data_to_coo_matrix_input_format(
        matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num
    )

    # Create the matrix from the sparse data representation generated above.
    print(f"Creating COO Matrix...")
    counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(features)), dtype=np.float32)

    # Convert the COO matrix to CSR for loading into AnnData
    print(f"Creating AnnData Matrix")
    
    adata = ad.AnnData(counts.tocsr(copy=False))
    adata.obs.index = original_cell_ids
    adata.var.index = feature_ids
    
    # See https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html?highlight=raw#anndata.AnnData.raw
    # for why we set 'raw' thusly: "The raw attribute is initialized with the current content of an object by setting:"
    adata.raw = adata
    
    print(f"Writing AnnData Matrix")    
    adata.write(Path(output), compression="gzip")
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

    counter = 0
    start = datetime.datetime.now()
    for row in matrix_data:
        cas_cell_index = row["cas_cell_index"]
        cas_feature_data = row["feature_data"]

        try:
            row_num = cas_cell_index_to_row_num[cas_cell_index]
        except KeyError as exc:
            raise Exception(
                f"ERROR: Unable to find entry for cas_cell_index: {cas_cell_index} in lookup table"
            ) from exc
        
        for e in cas_feature_data:
            cas_feature_index = e['feature_index']
            raw_counts = e['raw_counts']

            try:
                col_num = cas_feature_index_to_col_num[cas_feature_index]
            except KeyError as exc:
                raise Exception(
                    f"ERROR: Unable to find entry for cas_feature_index: {cas_feature_index} in lookup table"
                ) from exc
        
            rows.append(row_num)
            columns.append(col_num)
            data.append(raw_counts)
        
            counter = counter + 1
            if (counter % 1000000 == 0):
                end = datetime.datetime.now()
                elapsed = end - start
                print(f"Processed {counter} rows in {elapsed.total_seconds() * 1000}ms")
                start = end

    # TODO: do this with generators!!!
    return rows, columns, data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, description="Query CASP tables for random cells")
    parser.add_argument("--project", type=str, help="BigQuery Project", required=True)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=True)
    parser.add_argument("--extract_table_prefix", type=str, help="Prefix of extract tables", required=True)
    parser.add_argument("--start_bin", type=int, help="starting (inclusive) integer bin to extract", required=True)
    parser.add_argument("--end_bin", type=int, help="ending (inclusive) integer bin to extract", required=True)
    parser.add_argument("--output", type=str, help="Filename of the AnnData file", required=True)

    args = parser.parse_args()
    extract_minibatch_to_anndata(args.project, args.dataset, args.extract_table_prefix, args.start_bin, args.end_bin, args.output)
