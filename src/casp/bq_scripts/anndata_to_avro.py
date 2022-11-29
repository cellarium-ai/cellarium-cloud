"""
Converts input AnnData files to Avro format suitable for loading the Cell Annotation Service Pilot BigQuery schema
version 0.1.

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

"""

import argparse
import hashlib
import itertools
import json
import math
import multiprocessing
import os
import time

import anndata as ad
import numpy as np
from fastavro import parse_schema, writer
from google.api_core.exceptions import NotFound
from google.cloud import bigquery


def current_milli_time():
    """
    Return current time in millisecond precision.
    """
    return round(time.time() * 1000)


def write_avro(generator, parsed_schema, filename, progress_batch_size=10000):
    """
    Avro writing worker function. Calls back to `generator` for each row to be written.
    """
    start = current_milli_time()
    with open(filename, "a+b") as out:
        records = []
        counter = 0
        for record in generator():
            counter = counter + 1
            records.append(record)

            if counter % FLUSH_BATCH_SIZE == 0:
                writer(out, parsed_schema, records)
                records = []

            if counter % progress_batch_size == 0:
                end = current_milli_time()
                print(f"    Processed {counter} rows... in {end - start} ms")
                start = end

        if len(records) > 0:
            writer(out, parsed_schema, records)


def dump_core_matrix(row_array, col_array, data_array, row_offset, cas_cell_index_start, cas_feature_index_start, filename):
    """
    Write the raw X matrix.
    """
    schema = {
        "doc": "Raw data indexed by cell and feature id in the CAS BigQuery schema",
        "namespace": "cas",
        "name": "CellFeature",
        "type": "record",
        "fields": [
            {"name": "cas_cell_index", "type": "int"},
            {"name": "cas_feature_index", "type": "int"},
            {"name": "raw_counts", "type": "int"},
        ],
    }
    parsed_schema = parse_schema(schema)
    start = current_milli_time()
    #lg=row_lookup_d.get
    #rows = np.vectorize(lg)(row_array + row_offset)
    rows = row_array + row_offset + cas_cell_index_start
    print(f"    {filename} - Vectorize row lookup... in {current_milli_time() - start} ms")    
    start = current_milli_time()
    #cols = np.vectorize(col_lookup_d.get)(col_array)
    #cols = [ col_lookup[i] for i in col_array ]
    cols = col_array + cas_feature_index_start
    print(f"    {filename} - Processing {len(data_array)} elements")
    
    print(f"    Vectorize col lookup... in {current_milli_time() - start} ms")    
    start = current_milli_time()
    i_dat = data_array.astype(int)
    print(f"    {filename} - Convert types... in {current_milli_time() - start} ms")    

    # def raw_counts_generator():
    #     for row, col, dat in zip(rows, cols, i_dat):
    #         yield {
    #             "cas_cell_index": row,
    #             "cas_feature_index": col,
    #             "raw_counts": dat,
    #         }

    # write_avro(raw_counts_generator, parsed_schema, filename, progress_batch_size=1000000)
    # i =0
    # for row, col, dat in zip(rows, cols, i_dat):
    #     i = i + 1
    start = current_milli_time()
    with open(filename, "w") as out:
      for row, col, dat in zip(rows, cols, i_dat):
        out.write(f"{row},{col},{dat}\n")
    print(f"    {filename} - Write Chunk... in {current_milli_time() - start} ms")    

    # with open(filename, "a+b") as out:
    #      writer(out, parsed_schema, raw_counts_generator(), sync_interval=160000)


def dump_cell_info(adata, filename, cas_cell_index_start, ingest_id):
    """
    Write cell / obs data.
    """
    adata.obs["original_cell_id"] = adata.obs.index
    adata.obs["cas_cell_index"] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    # row_index_to_cas_cell_index = [None] * len(adata.obs)
    # for row in range(0, len(adata.obs)):
    #     row_index_to_cas_cell_index[row] = adata.obs["cas_cell_index"].iloc[[row]][0]

    schema = {
        "doc": "CAS Cell Info",
        "namespace": "cas",
        "name": "CellInfo",
        "type": "record",
        "fields": [
            {"name": "cas_cell_index", "type": "int"},
            {"name": "original_cell_id", "type": "string"},
            {"name": "cell_type", "type": "string"},
            {"name": "obs_metadata", "type": {"type": "string", "sqlType": "JSON"}},
            {"name": "cas_ingest_id", "type": "string"},
        ],
    }

    parsed_schema = parse_schema(schema)
    added_columns = ["cas_cell_index", "original_cell_id"]
    # A version of `obs` that does not contain our added entries, suitable for persisting to BigQuery.
    obs_original = adata.obs.drop(columns=added_columns)
    cells = adata.obs[added_columns + ["cell_type"]]

    def cell_generator():
        for idx, row in cells.iterrows():
            yield {
                "cas_cell_index": row["cas_cell_index"],
                "original_cell_id": row["original_cell_id"],
                "cell_type": row["cell_type"],
                "obs_metadata": obs_original.loc[idx].to_json(),
                "cas_ingest_id": ingest_id,
            }

    write_avro(cell_generator, parsed_schema, filename)

#    return row_index_to_cas_cell_index


def dump_feature_info(adata, filename, cas_feature_index_start, ingest_id):
    """
    Write feature / gene / var data.
    """
    adata.var["original_feature_id"] = adata.var.index
    adata.var["cas_feature_index"] = np.arange(cas_feature_index_start, cas_feature_index_start + len(adata.var))

    col_index_to_cas_feature_index = [None] * len(adata.var)
    for col in range(0, len(adata.var)):
        col_index_to_cas_feature_index[col] = adata.var["cas_feature_index"].iloc[[col]][0]

    schema = {
        "doc": "CAS Feature Info",
        "namespace": "cas",
        "name": "FeatureInfo",
        "type": "record",
        "fields": [
            {"name": "cas_feature_index", "type": "int"},
            {"name": "original_feature_id", "type": "string"},
            {"name": "feature_name", "type": "string"},
            {"name": "var_metadata", "type": {"type": "string", "sqlType": "JSON"}},
            {"name": "cas_ingest_id", "type": "string"},
        ],
    }

    parsed_schema = parse_schema(schema)
    added_columns = ["cas_feature_index", "original_feature_id"]
    # A version of `var` that does not contain our added entries, suitable for persisting to BigQuery.
    var_original = adata.var.drop(columns=added_columns)
    features = adata.var[added_columns + ["feature_name"]]

    def feature_generator():
        for i, row in features.iterrows():
            yield {
                "cas_feature_index": row["cas_feature_index"],
                "original_feature_id": row["original_feature_id"],
                "feature_name": row["feature_name"],
                "var_metadata": var_original.loc[i].to_json(),
                "cas_ingest_id": ingest_id,
            }

    write_avro(feature_generator, parsed_schema, filename)
    return col_index_to_cas_feature_index


def dump_ingest_info(adata, filename, ingest_id, load_uns_data):
    """
    Write ingest (AnnData-file level) data.
    """
    schema = {
        "doc": "CAS Ingest Info",
        "namespace": "cas",
        "name": "IngestInfo",
        "type": "record",
        "fields": [
            {"name": "cas_ingest_id", "type": "string"},
            {"name": "uns_metadata", "type": {"type": "string", "sqlType": "JSON"}},
            {"name": "ingest_timestamp", "type": ["null", "long"], "logicalType": ["null", "timestamp-millis"]},
        ],
    }

    parsed_schema = parse_schema(schema)

    # `uns` metadata is contained in a dict-like `OverloadedDict` type which does not offer a nice `tojson` method.
    # In practice there are often numpy value types embedded within `uns` which the `json` library does not handle
    # appropriately by default, hence the custom encoder class below.
    # https://stackoverflow.com/a/49677241
    class NumpyEncoder(json.JSONEncoder):
        """Special json encoder for numpy types"""

        def default(self, o):
            if isinstance(o, np.integer):
                return int(o)
            if isinstance(o, np.floating):
                return float(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            return json.JSONEncoder.default(self, o)

    def ingest_generator():
        # Some `uns` metadata has extremely large values that can break extract. Cap the allowed size of values and
        # substitute a `None` if values exceed the cap. This particular limit of 1 MiB was chosen somewhat arbitrarily;
        # it's plenty big but allows extracting the Tabula Sapiens endothlelial dataset that could not be extracted
        # without a cap.
        metadata_limit = 2**20
        uns = {}

        if load_uns_data is True:
            for idx, uncapped_val in adata.uns.data.items():
                uncapped_json = json.dumps(uncapped_val, cls=NumpyEncoder)
                if len(uncapped_json) > metadata_limit:
                    print(
                        f"AnnData `uns` contains a key `{idx}` whose JSONified value would have size {len(uncapped_json)} bytes."
                    )
                    print("Values this large can cause extraction to fail so this value is being nulled out.")
                    val = None
                else:
                    val = uncapped_val
                uns[idx] = val

        yield {"uns_metadata": json.dumps(uns, cls=NumpyEncoder), "cas_ingest_id": ingest_id, "ingest_timestamp": None}

    write_avro(ingest_generator, parsed_schema, filename)


def confirm_output_files_do_not_exist(filenames):
    """
    Throws if any of the specified files already exist.
    """
    existing = list(filter(os.path.exists, filenames))
    if len(existing) > 0:
        raise ValueError(
            f"Found existing output files, please rename or remove before running conversion: {', '.join(existing)}"
        )


# https://stackoverflow.com/a/3431838
def md5(filename):
    """
    Calculates the md5 hash of the specified file.
    """
    hash_md5 = hashlib.md5()
    chunks = 0
    progress_chunks = 2**18
    with open(filename, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            chunks += 1
            hash_md5.update(chunk)
            if chunks % progress_chunks == 0:
                print(f"    {chunks / progress_chunks} GiB...")
    return hash_md5.hexdigest()


def find_max_index(client, project, dataset, table, column):
    """
    Find the current maximum index in the specified table looking at the specified column.
    """
    dataset_id = f"{project}.{dataset}"
    table_id = f"{dataset_id}.{table}"
    try:
        _ = client.get_dataset(dataset_id)
    except NotFound as exc:
        raise ValueError(f"Dataset '{dataset_id}' not found, required to find max index in '{table_id}'.") from exc

    try:
        _ = client.get_table(table_id)
    except NotFound as exc:
        raise ValueError(f"Table '{table_id}' not found, required to find max index.") from exc

    query = f"""SELECT MAX({column}) AS max_id FROM `{dataset_id}.{table}`"""

    max_id = None
    job = client.query(query)
    for row in job.result():
        max_id = row["max_id"]

    # Default to -1 if no max id found. If the table is empty the query above will return a row with a null id.
    if max_id is None:
        max_id = -1

    return max_id


def process(input_file, cas_cell_index_start, cas_feature_index_start, avro_prefix, project, dataset, load_uns_data):
    """
    High level entry point, reads the input AnnData file and generates Avro files
    for ingest, cells, features, and raw / core data.
    """
    client = None
    if cas_cell_index_start is None:
        client = bigquery.Client(project=project)
        print("Looking for max id in `cas_cell_info`...")
        cas_cell_index_start = find_max_index(client, project, dataset, "cas_cell_info", "cas_cell_index") + 1
        print(f"cas_cell_index_start will be {cas_cell_index_start}")

    if cas_feature_index_start is None:
        if not client:
            client = bigquery.Client(project=project)
        print("Looking for max id in `cas_feature_info`...")
        cas_feature_index_start = find_max_index(client, project, dataset, "cas_feature_info", "cas_feature_index") + 1
        print(f"cas_feature_index_start will be {cas_feature_index_start}")

    avro_prefix = "cas" if not avro_prefix else avro_prefix
    print(f"Hashing input AnnData file '{input_file}' to generate ingest id...")
    ingest_id = f"cas-ingest-{md5(input_file)[:8]}"
    print(f"Generated ingest id '{ingest_id}'.")

    file_types = ["ingest_info", "cell_info", "feature_info", "raw_counts"]
    filenames = [f"{avro_prefix}_{file_type}.avro" for file_type in file_types]
    confirm_output_files_do_not_exist(filenames)
    (ingest_filename, cell_filename, feature_filename, raw_counts_filename) = filenames

    print(f"Loading input AnnData file '{input_file}'...")
    adata = ad.read(input_file, backed='r')

    print("Processing ingest metadata...")
    dump_ingest_info(adata, ingest_filename, ingest_id, load_uns_data)

    print("Processing cell/observation metadata...")
    dump_cell_info(adata, cell_filename, cas_cell_index_start, ingest_id)

    print("Processing feature/gene/variable metadata...")
    dump_feature_info(adata, feature_filename, cas_feature_index_start, ingest_id)

    print("Processing core data...")

    # recode the indexes to be the indexes of obs/var or should obs/var include these indices?
    # row_lookup_d = dict(enumerate(row_index_to_cas_cell_index))
    # col_lookup_d = dict(enumerate(col_index_to_cas_feature_index))

    total_cells = adata.raw.X.shape[0]
    batch_size = 10000
    num_batches = math.ceil(total_cells / batch_size)

    # ranges are start-inclusive and end-exclusive
    batches = [ (x, x*batch_size, min(total_cells, x*batch_size + batch_size))  for x in range(num_batches)]

    start = current_milli_time()
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        for (index, row_offset, end) in batches:
            chunk_raw_counts_filename = raw_counts_filename.replace(".avro",f".{index:06}.csv")
            args = [input_file, row_offset, end, cas_cell_index_start, cas_feature_index_start, chunk_raw_counts_filename]
            result = pool.apply_async(p_dump_core_matrix, args=(args,), error_callback=custom_error_callback)
            print(result)

        pool.close()
        pool.join()
    
    print(f"    Processed {total_cells} cells... in {current_milli_time() - start} ms")    

    print("Done.")

# handle any errors in the task function
def custom_error_callback(error):
    print(f'Got an Error: {error}', flush=True)

def p_dump_core_matrix(args):
    start = current_milli_time()

    (input_file, row_offset, end, cas_cell_index_start, cas_feature_index_start, raw_counts_filename) = args
    adata = ad.read(input_file, backed='r')
    coord = adata.raw.X[row_offset:end, :].tocoo()
    adata.file.close()
    print(f"    Read anndata, subset matrix... in {current_milli_time() - start} ms")    

    dump_core_matrix(coord.row, coord.col, coord.data, row_offset, cas_cell_index_start, cas_feature_index_start, raw_counts_filename)

# TODO: naive dump, compare
#    rows,cols = adata.X.nonzero()
#    for row,col in zip(rows,cols):
#        cas_cell_index = row_index_to_cas_cell_index[row]
#        cas_feature_index = col_index_to_cas_feature_index[col]
#        v = adata.X[row,col]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="Convert AnnData Single Cell Expression Data into format for loading into BQ"
    )

    parser.add_argument("--input", type=str, help="AnnData format input file", required=True)
    parser.add_argument(
        "--avro_prefix",
        type=str,
        help="Prefix to use for output Avro files, e.g. <avro_prefix>_cell_info.avro etc",
        required=False,
    )
    parser.add_argument("--cas_cell_index_start", type=int, help="starting number for cell index", required=False)
    parser.add_argument("--cas_feature_index_start", type=int, help="starting number for feature index", required=False)
    parser.add_argument("--flush_batch_size", type=int, help="max size of Avro batches to flush", required=False)
    parser.add_argument("--project", type=str, help="BigQuery Project", required=False)
    parser.add_argument("--dataset", type=str, help="BigQuery Dataset", required=False)
    parser.add_argument(
        "--load_uns_data", help="load uns (unstructured) metadata", default="False", action="store_true"
    )

    args = parser.parse_args()

    errors = []
    if args.cas_cell_index_start is None and not (args.project and args.dataset):
        errors.append(
            "if `--cas_cell_index_start` is not specified, --project and --dataset required to find start value"
        )

    if args.cas_feature_index_start is None and not (args.project and args.dataset):
        errors.append(
            "if `--cas_feature_index_start` is not specified, --project and --dataset required to find start value"
        )

    if errors:
        raise ValueError("\n\n".join(errors))

    FLUSH_BATCH_SIZE = 10000 if not args.flush_batch_size else args.flush_batch_size

    process(
        args.input,
        args.cas_cell_index_start,
        args.cas_feature_index_start,
        args.avro_prefix,
        args.project,
        args.dataset,
        args.load_uns_data,
    )
