"""
Selects a random subset of cells of a specified size from BigQuery and writes the data to output AnnData files.
"""

import concurrent.futures as concurrency
import datetime
import json
import multiprocessing
import os
import tempfile
import traceback
import typing as t
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import typer
from google.cloud import bigquery
from google.cloud.bigquery_storage import BigQueryReadClient, types
from google.oauth2.service_account import Credentials
from scipy.sparse import coo_matrix

from casp.bq_scripts import constants
from casp.bq_scripts.prepare_dataset_info import distinct_obs_column_values
from casp.data_manager import sql
from casp.services import utils


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


def add_cell_type_info(extract_bucket_path: str, bucket_name: str, adata: "ad.AnnData", cell_types: list):
    filename = "all_cell_types.csv"
    source_blob_name = f"{extract_bucket_path}/shared_meta/{filename}"

    if not os.path.exists(filename):
        utils.download_file_from_bucket(
            bucket_name=bucket_name, source_blob_name=source_blob_name, destination_file_name=filename
        )

    adata.obs["cell_type"] = cell_types
    df = pd.read_csv(filename, index_col=False)
    adata.obs.cell_type = pd.Categorical(adata.obs.cell_type.values, categories=df["cell_type"].values)


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

    adata.layers["measured_genes_mask"] = measured_genes_mask


def extract_minibatch_to_anndata(
    project: str,
    dataset: str,
    extract_table_prefix: str,
    start_bin: int,
    end_bin: int,
    output: str,
    bucket_name: str,
    extract_bucket_path: str,
    credentials: t.Optional[Credentials] = None,
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
    print("Extracting Feature Info...")

    features = get_features(project, dataset, f"{extract_table_prefix}__extract_feature_info", client)

    feature_ids = []
    cas_feature_index_to_col_num = {}
    for col_num, feature in enumerate(features):
        feature_ids.append(feature.original_feature_id)
        cas_feature_index_to_col_num[feature.cas_feature_index] = col_num

    print("Extracting Cell Info...")
    cells = get_cells_in_bin_range(
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=start_bin,
        end_bin=end_bin,
        client=client,
        obs_columns_to_include=obs_columns_to_include,
    )

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

    matrix_data = get_matrix_data(
        project, dataset, f"{extract_table_prefix}__extract_raw_count_matrix", start_bin, end_bin, credentials
    )

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
            extract_bucket_path=extract_bucket_path,
            bucket_name=bucket_name,
            adata=adata,
            cell_types=cell_types,
        )

    print("Adding Measured Genes Layer Mask...")
    add_measured_genes_layer_mask(
        extract_bucket_path=extract_bucket_path,
        bucket_name=bucket_name,
        adata=adata,
        cas_ingest_ids=cas_ingest_ids,
    )
    print("Adding Other Obs Columns...")
    for obs_column in obs_columns_values:
        adata.obs[obs_column] = obs_columns_values[obs_column]
        make_obs_categories_comprehensive(adata, obs_column, project, dataset, credentials)

    print("Writing AnnData Matrix")
    adata.write(Path(output), compression="gzip")
    print("Done.")


def make_obs_categories_comprehensive(adata: ad.AnnData, obs_column: str, project: str, dataset: str, credentials):
    dtype = adata.obs[obs_column].dtype
    if (dtype.name == "category") or (dtype.name == "object") or (dtype.name == "string"):
        print(f"Making '{obs_column}' ({dtype.name}) categories into a comprehensive categorical...")
        df = distinct_obs_column_values(obs_column=obs_column, project=project, dataset=dataset, credentials=credentials)
        adata.obs[obs_column] = pd.Categorical(adata.obs[obs_column].values, categories=df[obs_column].values)
    else:
        print(f"Skipping adata.obs['{obs_column}'] as it is a {dtype.name} type")


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


def extract_bin(
    project_id: str,
    dataset: str,
    extract_table_prefix: str,
    bin_number: int,
    file_name: str,
    file_path: str,
    output_bucket_name: str,
    extract_bucket_path: str,
    obs_columns_to_include: t.List[str],
) -> None:
    """
    Wrapper task `casp.bq_scripts.extract_minibatch_to_anndata` which processes exactly
    one bin at a time and saves the output anndata file to a GCS bucket.

    :param project_id: Google Cloud Project id
    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param bin_number: Bin to extract
    :param file_name: Name for a local file (.h5ad anndata) to save the output
    :param file_path: Local file path
    :param output_bucket_name: Name of GCS bucket
    :param extract_bucket_path: Path where the extract files and subdirectories should be located.
        Should correspond to the directory provided to prepare_extract script as current script uses `shared_meta`
        files.
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    try:
        extract_minibatch_to_anndata(
            project=project_id,
            dataset=dataset,
            extract_table_prefix=extract_table_prefix,
            start_bin=bin_number,
            end_bin=bin_number,
            output=file_path,
            bucket_name=output_bucket_name,
            extract_bucket_path=extract_bucket_path,
            obs_columns_to_include=obs_columns_to_include,
        )
        blob_name = f"{extract_bucket_path}/extract_files/{file_name}"
        utils.upload_file_to_bucket(local_file_name=file_path, blob_name=blob_name, bucket=output_bucket_name)
    except Exception as e:
        print(f"Error occurred in thread: {str(e)}")

    print(f"Processed bin {bin_number}")


def extract_bins_in_parallel_workers(
    project_id: str,
    dataset: str,
    extract_table_prefix: str,
    start_bin: int,
    end_bin: int,
    output_bucket_name: str,
    extract_bucket_path: str,
    obs_columns_to_include: str,
) -> None:
    """
    Wrapper task `casp.bq_scripts.extract_minibatch_to_anndata` which processes exactly
    one bin at a time and saves the output anndata file to a GCS bucket.

    :param project_id: Google Cloud Project id
    :param dataset: BigQuery Dataset
    :param extract_table_prefix: Prefix of extract tables
    :param start_bin: Start bin to extract
    :param end_bin: End bin to extract
    :param output_bucket_name: Name of GCS bucket
    :param extract_bucket_path: Path where the extract files and subdirectories should be located.
        Should correspond to the directory provided to prepare_extract script as current script uses `shared_meta`
        files.
    :param obs_columns_to_include: Optional list of columns from `cas_cell_info` table to include in ``adata.obs``.
        If not provided, no specific columns would be added to ``adata.obs`` apart from `cas_cell_index`.
        Note: It is required to provide the column names along with the aliases for the tables to which they belong.
        However, the output extract table would contain only the column names, without any aliases.
        Example: ``["c.cell_type", "c.donor_id", "c.sex", "i.dataset_id"]``
    """
    obs_columns_to_include = obs_columns_to_include.split(",")
    obs_columns_to_include = [x.split(".")[-1] for x in obs_columns_to_include]

    num_of_workers = multiprocessing.cpu_count()
    with tempfile.TemporaryDirectory() as temp_dir:
        with concurrency.ProcessPoolExecutor(max_workers=num_of_workers) as executor:
            futures = []
            for bin_number in range(start_bin, end_bin + 1, 1):
                file_name = f"extract_{bin_number}.h5ad"
                file_path = f"{temp_dir}/{file_name}"

                extract_task_kwargs = {
                    "project_id": project_id,
                    "dataset": dataset,
                    "extract_table_prefix": extract_table_prefix,
                    "bin_number": bin_number,
                    "file_name": file_name,
                    "file_path": file_path,
                    "output_bucket_name": output_bucket_name,
                    "extract_bucket_path": extract_bucket_path,
                    "obs_columns_to_include": obs_columns_to_include,
                }
                print("Submitting parallel extract job to thread executor...")
                future = executor.submit(extract_bin, **extract_task_kwargs)
                futures.append(future)

            done, not_done = concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)

            for future in done:
                try:
                    # Attempt to get the result of the future
                    _ = future.result()
                except Exception as e:
                    # If an exception is raised, print the exception details
                    print(f"Future: {future}")
                    print(f"Exception type: {type(e).__name__}")
                    print(f"Exception message: {e}")
                    # Format and print the full traceback
                    traceback.print_exception(type(e), e, e.__traceback__)
