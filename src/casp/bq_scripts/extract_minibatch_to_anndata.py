"""
Selects a random subset of cells of a specified size from BigQuery and writes the data to output AnnData files.
"""
import argparse
import datetime
import json
import os
import typing as t
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from google.cloud import bigquery
from google.cloud.bigquery_storage import BigQueryReadClient, types
from scipy.sparse import coo_matrix

from casp.services import utils

if t.TYPE_CHECKING:
    from google.oauth2.service_account import Credentials


class Feature:
    """
    Represents a CASP feature / gene.
    """

    def __init__(self, cas_feature_index, original_feature_id):
        self.cas_feature_index = cas_feature_index
        self.original_feature_id = original_feature_id


def get_features(project, dataset, extract_feature_table, client):
    """
    Retrieve a list of all feature objects from the specified extract feature table, ordered by index from schema table.
    """
    sql = f"""
    SELECT cas_feature_index, original_feature_id FROM
        `{project}.{dataset}.{extract_feature_table}`
    ORDER BY index
    """

    result = client.query(sql)
    features = [Feature(row["cas_feature_index"], row["original_feature_id"]) for row in result]
    return features


def get_cells_in_bin_range(
    project: str, dataset: str, extract_table_prefix: str, start_bin: int, end_bin: int, client: "bigquery.Client"
) -> t.Iterable[t.Dict[str, t.Any]]:
    """
    Retrieve a list of all cells between the specified start and end bins (both inclusive), ordered by cas_cell_index.
    """
    sql = f"""
    SELECT cas_cell_index,
           cas_ingest_id,
           cell_type,
           total_mrna_umis
    FROM `{project}.{dataset}.{extract_table_prefix}__extract_cell_info`
    WHERE extract_bin BETWEEN {start_bin} AND {end_bin}
    """

    return [x for x in client.query(sql)]


def get_matrix_data(project, dataset, table, start_bin, end_bin, credentials=None):
    if credentials is None:
        client = BigQueryReadClient()
    else:
        client = BigQueryReadClient(credentials=credentials)

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


def add_cell_type_info(
    dataset: str, extract_table_prefix: str, bucket_name: str, adata: "ad.AnnData", cell_types: list
):
    filename = "all_cell_types.csv"
    source_blob_name = f"{dataset}_{extract_table_prefix}_info/{filename}"

    if not os.path.exists(filename):
        utils.download_file_from_bucket(
            bucket_name=bucket_name, source_blob_name=source_blob_name, destination_file_name=filename
        )

    adata.obs["cell_type"] = cell_types
    df = pd.read_csv(filename, index_col=False)
    adata.obs.cell_type = pd.Categorical(adata.obs.cell_type.values, categories=df["cell_type"].values)


def add_measured_genes_layer_mask(
    dataset: str, extract_table_prefix: str, bucket_name: str, adata: "ad.AnnData", cas_ingest_ids: list
):
    filename = "measured_genes_info.csv"
    source_blob_name = f"{dataset}_{extract_table_prefix}_info/{filename}"

    if not os.path.exists(filename):
        utils.download_file_from_bucket(
            bucket_name=bucket_name, source_blob_name=source_blob_name, destination_file_name=filename
        )

    adata.obs["cas_ingest_id"] = cas_ingest_ids
    df = pd.read_csv(filename, index_col="cas_ingest_id").astype(bool)

    measured_genes_mask = df.loc[adata.obs["cas_ingest_id"], adata.var.index].to_numpy()

    adata.layers["measured_genes_mask"] = measured_genes_mask


def extract_minibatch_to_anndata(
    project,
    dataset,
    extract_table_prefix,
    start_bin,
    end_bin,
    output,
    bucket_name,
    credentials: "Credentials" = None,
    obs_columns_to_include: t.Optional[t.List[str]] = None,
):
    """
    Main function to extract a minibatch from a prepared training extract and write the associated data
    to an AnnData file.

    :param project: BigQuery Project
    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param start_bin: Starting (inclusive) integer bin to extract
    :param end_bin: Ending (inclusive) integer bin to extract
    :param output: Filename of the AnnData file
    :param bucket_name: Bucket name where to save extracted minibatch
    :param credentials: Google Cloud Credentials with the access to BigQuery
    :param obs_columns_to_include: Columns to include in `obs` data frame in extracted `anndata.AnnData` file
        Mapped from extract query output
    """
    if credentials is None:
        client = bigquery.Client(project=project)
    else:
        client = bigquery.Client(project=project, credentials=credentials)

    print(f"Getting extract bin {start_bin}-{end_bin} data from {project}.{dataset}.{extract_table_prefix}*...")

    # Read the feature information and store for later.
    print("Extracting Feature Info...")

    features = get_features(project, dataset, f"{extract_table_prefix}__extract_feature_info", client)

    feature_ids = []
    cas_feature_index_to_col_num = {}
    for col_num, feature in enumerate(features):
        feature_ids.append(feature.original_feature_id)
        cas_feature_index_to_col_num[feature.cas_feature_index] = col_num

    print("Extracting Cell Info...")
    cells = get_cells_in_bin_range(project, dataset, extract_table_prefix, start_bin, end_bin, client)

    original_cell_ids = []
    cas_ingest_ids = []
    obs_columns_values = {k: [] for k in obs_columns_to_include}
    cas_cell_index_to_row_num = {}
    for row_num, cell in enumerate(cells):
        original_cell_ids.append(cell["cas_cell_index"])
        cas_ingest_ids.append(cell["cas_ingest_id"])
        cas_cell_index_to_row_num[cell["cas_cell_index"]] = row_num
        # Mapping additional obs columns
        for obs_column in obs_columns_to_include:
            obs_columns_values[obs_column].append(cell[obs_column])

    matrix_data = get_matrix_data(
        project, dataset, f"{extract_table_prefix}__extract_raw_count_matrix", start_bin, end_bin, credentials
    )
    # print(next(iter(matrix_data)))
    print("Converting Matrix Data to COO format...")

    rows, columns, data = convert_matrix_data_to_coo_matrix_input_format(
        matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num
    )

    # Create the matrix from the sparse data representation generated above.
    print("Creating COO Matrix...")
    counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(features)), dtype=np.float32)

    # Convert the COO matrix to CSR for loading into AnnData
    print("Creating AnnData Matrix")

    adata = ad.AnnData(counts.tocsr(copy=False))
    adata.obs.index = original_cell_ids
    adata.var.index = feature_ids

    print("Adding Cell Types Categorical Info...")
    if "cell_type" in obs_columns_values:
        # To avoid adding and overwriting cell types, pop them from dictionary:
        cell_types = obs_columns_values.pop("cell_type")
        add_cell_type_info(
            dataset=dataset,
            extract_table_prefix=extract_table_prefix,
            bucket_name=bucket_name,
            adata=adata,
            cell_types=cell_types,
        )

    print("Adding Measured Genes Layer Mask...")
    add_measured_genes_layer_mask(
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        bucket_name=bucket_name,
        adata=adata,
        cas_ingest_ids=cas_ingest_ids,
    )
    print("Adding Other Obs Columns...")
    for obs_column in obs_columns_values:
        adata.obs[obs_column] = obs_columns_values[obs_column]

    print("Writing AnnData Matrix")
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
            cas_feature_index = e["feature_index"]
            raw_counts = e["raw_counts"]

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
            if counter % 1000000 == 0:
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
    extract_minibatch_to_anndata(
        args.project, args.dataset, args.extract_table_prefix, args.start_bin, args.end_bin, args.output
    )
