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
import scipy.sparse as sp
from google.cloud import bigquery
from google.cloud.bigquery_storage import BigQueryReadClient, types
from scipy.sparse import coo_matrix

from casp.bq_scripts import constants
from casp.bq_scripts.prepare_dataset_info import distinct_obs_column_values
from casp.data_manager import sql
from casp.services import utils

if t.TYPE_CHECKING:
    from google.oauth2.service_account import Credentials


def get_cells_in_bin_range(
    project: str,
    dataset: str,
    extract_table_prefix: str,
    start_bin: int,
    end_bin: int,
    client: "bigquery.Client",
    obs_columns_to_include: t.Optional[t.List[str]] = None,
) -> t.Iterable[t.Dict[str, t.Any]]:
    """
    Retrieve a list of all cells between the specified start and end bins (both inclusive), ordered by cas_cell_index.
    """
    template_data = sql.TemplateData(
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        select=obs_columns_to_include,
        start_bin=start_bin,
        end_bin=end_bin,
    )
    rendered_sql = sql.render(template_path=constants.GET_CELLS_IN_BIN_RANGE_SQL_DIR, template_data=template_data)
    return [x for x in client.query(rendered_sql)]


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
    requested_session.data_format = types.DataFormat.ARROW

    requested_session.read_options.selected_fields = ["cas_cell_index", "cas_feature_index", "raw_counts"]
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
    # reader = client.read_rows(session.streams[0].name)

    # rows = reader.rows(session)
    # return rows
    return session


def get_var_dataframe(project, dataset, extract_feature_table, client) -> pd.DataFrame:
    """
    Retrieve a dataframe of all feature objects from the specified extract feature table, ordered by index from schema table.
    """
    sql = f"""
        WITH ranked_feature_info AS (
            SELECT
                original_feature_id,
                feature_name,
                feature_biotype,
                feature_reference,
                ROW_NUMBER() OVER (PARTITION BY original_feature_id ORDER BY feature_name) AS row_num
            FROM
                `{project}.{dataset}.cas_feature_info`
        )
        SELECT 
            extracted_features.cas_feature_index, 
            extracted_features.original_feature_id, 
            ranked_feature_info.feature_name, 
            ranked_feature_info.feature_biotype, 
            ranked_feature_info.feature_reference 
        FROM
            `{project}.{dataset}.{extract_feature_table}` extracted_features
        INNER JOIN
            ranked_feature_info
        ON
            extracted_features.original_feature_id = ranked_feature_info.original_feature_id
        WHERE
        ranked_feature_info.row_num = 1
        ORDER BY cas_feature_index
    """
    return client.query(sql).to_dataframe()


def get_feature_lookup(project: str, dataset: str, var: pd.DataFrame, client):
    """
    Retrieve a dictionary mapping cas_raw_count_matrix.cas_feature_index 
    to {extract_table}__extract_feature_info.cas_feature_index
    thereby enabling us to collapse the same features based on our schema.
    """
    # sql = f"""
    #     WITH feature_definitions AS (
    #         SELECT
    #             cas_feature_index,
    #             original_feature_id
    #         FROM
    #             `{project}.{dataset}.cas_feature_info`
    #     ), 
    #     feature_extract AS (
    #         SELECT
    #             cas_feature_index,
    #             original_feature_id
    #         FROM
    #             `{project}.{dataset}.{extract_feature_table}`
    #     )
    #     SELECT
    #         count_matrix.cas_feature_index AS count_matrix_cas_feature_index,
    #         feature_definitions.original_feature_id,
    #         feature_extract.cas_feature_index AS extract_cas_feature_index
    #     FROM
    #         `{project}.{dataset}.cas_raw_count_matrix` 
    #     AS count_matrix
    #     INNER JOIN 
    #         feature_definitions
    #     ON
    #         count_matrix.cas_feature_index = feature_definitions.cas_feature_index
    #     INNER JOIN
    #         feature_extract
    #     ON
    #         feature_definitions.original_feature_id = feature_extract.original_feature_id
    # """
    sql = f"""
        SELECT
            cas_feature_index,
            original_feature_id
        FROM
            `{project}.{dataset}.cas_feature_info`
    """
    df = client.query(sql).to_dataframe()
    df_var = pd.DataFrame(data={'original_feature_id': var.index.values, 'var_index': range(len(var))})
    df = pd.merge(left=df, right=df_var, on='original_feature_id', how='left')
    lookup_dict = (
        df[['cas_feature_index', 'var_index']]
        .set_index('cas_feature_index')
        ['var_index']
        .to_dict()
    )
    return lookup_dict


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


def make_obs_categories_comprehensive(adata: ad.AnnData, obs_column: str, project: str, dataset: str, credentials):
    dtype = adata.obs[obs_column].dtype
    if (dtype.name == "category") or (dtype.name == "object") or (dtype.name == "string"):
        print(f"Making '{obs_column}' ({dtype.name}) categories into a comprehensive categorical...")
        try:
            df = distinct_obs_column_values(obs_column=obs_column, project=project, dataset=dataset, credentials=credentials)
            adata.obs[obs_column] = pd.Categorical(adata.obs[obs_column].values, categories=df[obs_column].values)
        except Exception as e:
            print(f"WARNING!: Failed to make '{obs_column}' ({dtype.name}) categories comprehensive:\n{str(e)}")
    else:
        print(f"Skipping adata.obs['{obs_column}'] as it is a {dtype.name} type")


def add_measured_genes_layer_mask(extract_bucket_path, bucket_name: str, adata: "ad.AnnData", cas_ingest_ids: list):
    filename = "measured_genes_info.csv"
    source_blob_name = f"{extract_bucket_path}/shared_meta/{filename}"

    if not os.path.exists(filename):
        utils.download_file_from_bucket(
            bucket_name=bucket_name, source_blob_name=source_blob_name, destination_file_name=filename
        )

    adata.obs["cas_ingest_id"] = cas_ingest_ids
    df = pd.read_csv(filename, index_col="cas_ingest_id").astype(bool)

    measured_genes_mask = df.loc[adata.obs["cas_ingest_id"], adata.var.index].to_numpy()

    adata.layers["feature_not_measured"] = sp.csr_matrix(np.logical_not(measured_genes_mask))


def extract_minibatch_to_anndata(
    project: str,
    dataset: str,
    extract_table_prefix: str,
    start_bin: int,
    end_bin: int,
    bucket_name: str,
    extract_bucket_path: str,
    credentials: "Credentials" = None,
    obs_columns_to_include: t.Optional[t.List[str]] = None,
) -> ad.AnnData:
    """
    Main function to extract a minibatch from a prepared training extract and write the associated data
    to an AnnData file.

    :param project: BigQuery Project
    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param start_bin: Starting (inclusive) integer bin to extract
    :param end_bin: Ending (inclusive) integer bin to extract
    :param bucket_name: Bucket name where to save extracted minibatch
    :param extract_bucket_path: Path where the extract files and subdirectories should be located.
        Should correspond to the directory provided to prepare_extract script as current script uses `shared_meta`
        files.
    :param credentials: Google Cloud Credentials with the access to BigQuery
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    if credentials is None:
        client = bigquery.Client(project=project)
    else:
        client = bigquery.Client(project=project, credentials=credentials)

    print(f"Getting extract bin {start_bin}-{end_bin} data from {project}.{dataset}.{extract_table_prefix}*...")

    # Read the feature information and store for later.
    # TODO: we do not need to do this every time, do we?  this happens over and over and over.
    # NOTE: so I did try to factor this out but my implementation did not work well, presumably because 
    # I was passing too much data around to subprocesses
    # NOTE: copilot thinks this is a good candidate for a cache
    print(f"bin {start_bin} Extracting Feature Info...")
    feature_df = get_var_dataframe(project, dataset, f"{extract_table_prefix}__extract_feature_info", client)
    var = feature_df.set_index("original_feature_id").drop(columns=["cas_feature_index"])
    cas_feature_index_to_col_num = get_feature_lookup(project, dataset, var, client)
    # cas_feature_index_to_col_num = dict(zip(feature_df["cas_feature_index"], range(len(feature_df))))

    print(f"bin {start_bin} Extracting Cell Info...")
    cells = get_cells_in_bin_range(
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=start_bin,
        end_bin=end_bin,
        client=client,
        obs_columns_to_include=obs_columns_to_include,
    )
    print(f"... bin {start_bin} grabbed {len(cells)} cells from extract_cell_info table")

    # Getting rid of table aliases in `obs_columns_to_include`
    _obs_columns_to_include = sql.mako_helpers.remove_leading_alias(column_names=obs_columns_to_include)
    original_cell_ids = []
    cas_ingest_ids = []
    obs_columns_values = {k: [] for k in _obs_columns_to_include}
    cas_cell_index_to_row_num = {}
    for row_num, cell in enumerate(cells):
        original_cell_ids.append(cell["cas_cell_index"])
        cas_ingest_ids.append(cell["cas_ingest_id"])
        cas_cell_index_to_row_num[cell["cas_cell_index"]] = row_num
        # Mapping additional obs columns
        for obs_column in _obs_columns_to_include:
            obs_columns_values[obs_column].append(cell[obs_column])

    print(f"bin {start_bin} Extracting Count Matrix...")
    matrix_data = get_matrix_data(
        project, dataset, f"{extract_table_prefix}__extract_raw_count_matrix", start_bin, end_bin, credentials
    )

    print(f"bin {start_bin} Converting Matrix Data to COO format...")

    rows, columns, data = convert_matrix_data_to_coo_matrix_input_format(
        matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num, credentials=credentials,
    )

    print(f'bin {start_bin} max value in rows is {max(rows)}')
    print(f'bin {start_bin} max value in columns is {max(columns)}')
    print(f'bin {start_bin} max value in data is {max(data)}')
    print(f'bin {start_bin} shape is {(len(cells), len(var))}')

    # Create the matrix from the sparse data representation generated above.
    print(f"bin {start_bin} Creating COO Matrix...")
    counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(var)), dtype=np.float32)

    # Convert the COO matrix to CSR for loading into AnnData
    print(f"bin {start_bin} Creating AnnData Matrix")

    adata = ad.AnnData(X=counts.tocsr(copy=False), var=var)
    adata.obs.index = original_cell_ids
    adata.obs.index.name = 'cas_cell_index'

    print(f"bin {start_bin} Adding Measured Genes Layer Mask...")
    add_measured_genes_layer_mask(
        extract_bucket_path=extract_bucket_path,
        bucket_name=bucket_name,
        adata=adata,
        cas_ingest_ids=cas_ingest_ids,
    )

    print(f"bin {start_bin} Adding Other Obs Columns...")
    for obs_column in obs_columns_values:
        adata.obs[obs_column] = obs_columns_values[obs_column]
        make_obs_categories_comprehensive(adata, obs_column, project, dataset, credentials)

    return adata

    # print("Writing AnnData Matrix")
    # adata.write(Path(output), compression="gzip")
    # print("Done.")

    # return Path(output)


def convert_matrix_data_to_coo_matrix_input_format(
    matrix_data, cas_cell_index_to_row_num, cas_feature_index_to_col_num, credentials,
):
    """
    Convert raw matrix data to coordinate format suitable for writing out an AnnData file.
    """
    # rows = []
    # columns = []
    # data = []

    # # TODO: can you vectorize this thing using pandas?  make a dataframe and then do 
    # # row['cas_cell_index'].map(cas_cell_index_to_row_num) to get the row numbers ?  let pandas use C backend

    # counter = 0
    start = datetime.datetime.now()
    # for row in matrix_data:
    #     cas_cell_index = row["cas_cell_index"]
    #     cas_feature_data = row["feature_data"]

    #     try:
    #         row_num = cas_cell_index_to_row_num[cas_cell_index]
    #     except KeyError as exc:
    #         raise Exception(
    #             f"ERROR: Unable to find entry for cas_cell_index: {cas_cell_index} in lookup table"
    #         ) from exc

    #     for e in cas_feature_data:
    #         cas_feature_index = e["feature_index"]
    #         raw_counts = e["raw_counts"]

    #         try:
    #             col_num = cas_feature_index_to_col_num[cas_feature_index]
    #         except KeyError as exc:
    #             raise Exception(
    #                 f"ERROR: Unable to find entry for cas_feature_index: {cas_feature_index} in lookup table"
    #             ) from exc

    #         rows.append(row_num)
    #         columns.append(col_num)
    #         data.append(raw_counts)

    #         counter = counter + 1
    #         if counter % 1000000 == 0:
    #             end = datetime.datetime.now()
    #             elapsed = end - start
    #             print(f"Processed {counter} rows in {elapsed.total_seconds() * 1000}ms")
    #             start = end

    # # TODO: do this with generators!!!
    # return rows, columns, data
    if credentials is None:
        client = BigQueryReadClient()
    else:
        client = BigQueryReadClient(credentials=credentials)

    dfs = []
    for i, stream in enumerate(matrix_data.streams):
        reader = client.read_rows(stream.name)
        for batch in reader.rows(matrix_data).pages:
            arrow_table = batch.to_arrow()
            df_batch = arrow_table.to_pandas()
            dfs.append(df_batch)

        end = datetime.datetime.now()
        elapsed = end - start
        print(f"Processed stream {i} in {elapsed.total_seconds()} sec")
        start = end

    df = pd.concat(dfs, ignore_index=True)
    df["row"] = df["cas_cell_index"].map(cas_cell_index_to_row_num)
    df["col"] = df["cas_feature_index"].map(cas_feature_index_to_col_num)
    return df["row"].values, df["col"].values, df["raw_counts"].values


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
