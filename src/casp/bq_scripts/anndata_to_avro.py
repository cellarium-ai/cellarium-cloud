"""
Converts input AnnData files to Avro format suitable for loading the Cell Annotation Service Pilot BigQuery schema
version 0.1.

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

"""

import gc
import hashlib
import json
import math
import multiprocessing
import os
import time
import typing as t

import h5py
import numpy as np
from anndata._core.anndata import AnnData
from anndata._io.h5ad import _clean_uns, _read_raw
from anndata._io.specs import read_elem
from fastavro import parse_schema, writer
from google.api_core.exceptions import NotFound
from google.cloud import bigquery
from scipy import sparse


def current_milli_time():
    """
    Return current time in millisecond precision.
    """
    return round(time.time() * 1000)


# Default value for maximum batch size of avro files when they are being created by avro writer
FLUSH_BATCH_SIZE_DEFAULT = 10000
# Default value for count matrix multiprocessing batch size
COUNT_MATRIX_MULTIPROCESSING_BATCH_SIZE_DEFAULT = 5000
# Default value for feature id lookup.
ORIGINAL_FEATURE_ID_LOOKUP_DEFAULT = "index"


def write_avro(
    generator, parsed_schema, filename, progress_batch_size=10000, flush_batch_size: int = FLUSH_BATCH_SIZE_DEFAULT
):
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

            if counter % flush_batch_size == 0:
                writer(out, parsed_schema, records)
                records = []

            if counter % progress_batch_size == 0:
                end = current_milli_time()
                print(f"    Processed {counter} rows... in {end - start} ms")
                start = end

        if len(records) > 0:
            writer(out, parsed_schema, records)


def get_x_matrix(adata: AnnData) -> t.Optional[t.Union[np.ndarray, sparse.spmatrix]]:
    """
    If there `raw` is not None use `raw.X` otherwise use just X
    :param adata: Anndata file to use
    """
    if adata.raw is not None:
        return adata.raw.X

    return adata.X


def dump_core_matrix(
    row_array, col_array, data_array, row_offset, cas_cell_index_start, cas_feature_index_start, filename
):
    """
    Write the raw X matrix.
    """
    start = current_milli_time()
    # lg=row_lookup_d.get
    # rows = np.vectorize(lg)(row_array + row_offset)
    rows = row_array + row_offset + cas_cell_index_start
    print(f"    {filename} - Vectorize row lookup... in {current_milli_time() - start} ms")
    start = current_milli_time()
    # cols = np.vectorize(col_lookup_d.get)(col_array)
    # cols = [ col_lookup[i] for i in col_array ]
    cols = col_array + cas_feature_index_start
    print(f"    {filename} - Processing {len(data_array)} elements")

    print(f"    {filename} - Vectorize col lookup... in {current_milli_time() - start} ms")
    start = current_milli_time()
    i_dat = data_array.astype(int)
    print(f"    {filename} - Convert types... in {current_milli_time() - start} ms")

    start = current_milli_time()
    with open(filename, "w") as out:
        for row, col, dat in zip(rows, cols, i_dat):
            out.write(f"{row},{col},{dat}\n")
    print(f"    {filename} - Write Chunk... in {current_milli_time() - start} ms")


def dump_cell_info(adata: AnnData, filename: str, cas_cell_index_start: int, ingest_id: str):
    """
    Write cell / obs data.
    """
    adata.obs["original_cell_id"] = adata.obs.index
    adata.obs["cas_cell_index"] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    schema = {
        "doc": "CAS Cell Info",
        "namespace": "cas",
        "name": "CellInfo",
        "type": "record",
        "fields": [
            {"name": "cas_cell_index", "type": "int"},
            {"name": "original_cell_id", "type": "string"},
            {"name": "assay_ontology_term_id", "type": "string"},
            {"name": "cell_type_ontology_term_id", "type": "string"},
            {"name": "development_stage_ontology_term_id", "type": "string"},
            {"name": "disease_ontology_term_id", "type": "string"},
            {"name": "donor_id", "type": "string"},
            {"name": "is_primary_data", "type": "boolean"},
            {"name": "organism_ontology_term_id", "type": "string"},
            {"name": "self_reported_ethnicity_ontology_term_id", "type": "string"},
            {"name": "sex_ontology_term_id", "type": "string"},
            {"name": "suspension_type", "type": "string"},
            {"name": "tissue_ontology_term_id", "type": "string"},
            {"name": "assay", "type": "string"},
            {"name": "cell_type", "type": "string"},
            {"name": "development_stage", "type": "string"},
            {"name": "disease", "type": "string"},
            {"name": "organism", "type": "string"},
            {"name": "self_reported_ethnicity", "type": "string"},
            {"name": "sex", "type": "string"},
            {"name": "tissue", "type": "string"},
            {"name": "obs_metadata_extra", "type": {"type": "string", "sqlType": "JSON"}},
            {"name": "cas_ingest_id", "type": "string"},
            {"name": "total_mrna_umis", "type": ["null", "int"], "default": None},
        ],
    }

    parsed_schema = parse_schema(schema)
    added_columns = [
        "cas_cell_index",
        "original_cell_id",
        "assay_ontology_term_id",
        "cell_type_ontology_term_id",
        "development_stage_ontology_term_id",
        "disease_ontology_term_id",
        "donor_id",
        "is_primary_data",
        "organism_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "sex_ontology_term_id",
        "suspension_type",
        "tissue_ontology_term_id",
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "organism",
        "self_reported_ethnicity",
        "sex",
        "tissue",
    ]
    # A version of `obs` that does not contain our added entries, suitable for persisting to BigQuery.
    obs_original = adata.obs.drop(columns=added_columns)
    cells = adata.obs[added_columns]

    def cell_generator():
        for idx, row in cells.iterrows():
            obs_data_to_include = {}
            for column_to_add in added_columns:
                obs_data_to_include[column_to_add] = row.get(column_to_add)

            yield {
                **obs_data_to_include,
                "obs_metadata_extra": obs_original.loc[idx].dropna().to_json(),
                "cas_ingest_id": ingest_id,
            }

    write_avro(cell_generator, parsed_schema, filename)


def dump_feature_info(
    adata, filename, cas_feature_index_start, ingest_id, original_feature_id_lookup=ORIGINAL_FEATURE_ID_LOOKUP_DEFAULT
):
    """
    Write feature / gene / var data.
    """
    if original_feature_id_lookup == ORIGINAL_FEATURE_ID_LOOKUP_DEFAULT:
        adata.var["original_feature_id"] = adata.var.index
    else:
        adata.var["original_feature_id"] = adata.var[original_feature_id_lookup]

    adata.var["cas_feature_index"] = np.arange(cas_feature_index_start, cas_feature_index_start + len(adata.var))

    schema = {
        "doc": "CAS Feature Info",
        "namespace": "cas",
        "name": "FeatureInfo",
        "type": "record",
        "fields": [
            {"name": "cas_feature_index", "type": "int"},
            {"name": "original_feature_id", "type": "string"},
            {"name": "feature_name", "type": "string"},
            {"name": "feature_biotype", "type": "string"},
            {"name": "feature_is_filtered", "type": "boolean"},
            {"name": "feature_reference", "type": "string"},
            {"name": "var_metadata_extra", "type": {"type": "string", "sqlType": "JSON"}},
            {"name": "cas_ingest_id", "type": "string"},
        ],
    }

    parsed_schema = parse_schema(schema)
    added_columns = [
        "cas_feature_index",
        "original_feature_id",
        "feature_name",
        "feature_biotype",
        "feature_is_filtered",
        "feature_reference",
    ]
    # A version of `var` that does not contain our added entries, suitable for persisting to BigQuery.
    var_original = adata.var.drop(columns=added_columns)
    features = adata.var[added_columns]

    def feature_generator():
        for i, row in features.iterrows():
            yield {
                "cas_feature_index": row["cas_feature_index"],
                "original_feature_id": row["original_feature_id"],
                # Try to retrieve these attributes as well, if not presented put None
                "feature_name": row.get("feature_name"),
                "feature_biotype": row.get("feature_biotype"),
                "feature_is_filtered": row.get("feature_is_filtered"),
                "feature_reference": row.get("feature_reference"),
                "var_metadata_extra": var_original.loc[i].dropna().to_json(),
                "cas_ingest_id": ingest_id,
            }

    write_avro(feature_generator, parsed_schema, filename)
    # return col_index_to_cas_feature_index


def dump_ingest_info(
    adata: AnnData,
    filename: str,
    dataset_id: str,
    dataset_version_id: str,
    source_filename: str,
    ingest_id: str,
    load_uns_data: bool,
    prefix: str,
    included_adata_uns_keys: t.List[str] = None,
):
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
            {"name": "dataset_id", "type": "string"},
            {"name": "dataset_version_id", "type": "string"},
            {"name": "source_filename", "type": "string"},
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
            else:
                try:
                    return json.JSONEncoder.default(self, o)
                except TypeError:
                    # In case if the instance is still not serializable just keep track of its type
                    return o.__class__.__name__

    def ingest_generator():
        # Some `uns` metadata has extremely large values that can break extract. Cap the allowed size of values and
        # substitute a `None` if values exceed the cap. This particular limit of 1 MiB was chosen somewhat arbitrarily;
        # it's plenty big but allows extracting the Tabula Sapiens endothlelial dataset that could not be extracted
        # without a cap.
        metadata_limit = 2**20
        uns = {}

        if load_uns_data:
            for idx, uncapped_val in adata.uns.data.items():
                if included_adata_uns_keys is not None:
                    if idx not in included_adata_uns_keys:
                        continue

                uncapped_json = json.dumps(uncapped_val, cls=NumpyEncoder)
                if len(uncapped_json) > metadata_limit:
                    print(
                        f"AnnData `uns` contains a key `{idx}` whose JSONified value would have size {len(uncapped_json)} bytes."
                    )
                    print("Values this large can cause extraction to fail so this value is being nulled out.")
                else:
                    pass
                val = uncapped_val
                uns[idx] = val

        yield {
            "uns_metadata": json.dumps(uns, cls=NumpyEncoder),
            "dataset_id": dataset_id,
            "dataset_version_id": dataset_version_id,
            "source_filename": source_filename,
            "cas_ingest_id": ingest_id,
            "ingest_timestamp": None,
            "ingest_prefix_name": prefix,
        }

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


def process(
    input_file: str,
    prefix: str,
    czi_dataset_id: str,
    czi_dataset_version_id: str,
    project: str,
    load_uns_data: bool,
    dataset: t.Optional[str] = None,
    cas_cell_index_start: t.Optional[int] = None,
    cas_feature_index_start: t.Optional[int] = None,
    included_adata_uns_keys: t.List = None,
    original_feature_id_lookup=ORIGINAL_FEATURE_ID_LOOKUP_DEFAULT,
    count_matrix_multiprocessing_batch_size: int = COUNT_MATRIX_MULTIPROCESSING_BATCH_SIZE_DEFAULT,
    save_directory: t.Optional[str] = None,
):
    """
    High level entry point, reads the input AnnData file and generates Avro files
    for ingest, cells, features, and raw / core data.

    :param original_feature_id_lookup: A column name in var dataframe from where to get original feature ids.
    In most of the cases it will be a column with ENSEMBL gene IDs. Default is `index` which means that
    an index column of var dataframe would be used.

    :param input_file: Input file that has to be ingested
    :param prefix: Prefix for filenames
    :param czi_dataset_id: Dataset id from cellxgene census
    :param czi_dataset_version_id: Version id of the current dataset from cellxgene census
    :param project: Google Cloud project which is used for connecting to BigQuery
    :param load_uns_data: Whether to load uns (unstructured) metadata
    :param dataset: BigQuery dataset. Used when max indexes for cells and indexes has to be retrieved from BigQuery
    :param cas_cell_index_start: Starting number for cell index. If ``None``, increment of maximum index would be used.
        ``None`` requires ``dataset`` to be set as it will use dataset to get the maximum index in the cell_info table
    :param cas_feature_index_start: Starting number for feature index. If ``None``, increment of maximum index would be
        used. ``None`` requires ``dataset`` to be set as it will use dataset to get the maximum index in the cell_info
        table
    :param included_adata_uns_keys: List with a set of keys that need to be dumped in ingest. If None, dump all.
    :param count_matrix_multiprocessing_batch_size: Size of count matrix csv files batches that are used in ingest
    :param save_directory: Directory where to save local. If ``None`` current local directory will be used.
        It is recommended to use a temporary folder for this processing.
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

    prefix = "cas" if not prefix else prefix
    print(f"Hashing input AnnData file '{input_file}' to generate ingest id...")
    ingest_id = f"cas-ingest-{md5(input_file)[:8]}"
    print(f"Generated ingest id '{ingest_id}'.")

    file_types = ["ingest_info", "cell_info", "feature_info", "raw_counts"]
    filenames = [f"{prefix}_{file_type}.avro" for file_type in file_types]
    confirm_output_files_do_not_exist(filenames)

    directory_prefix = f"{save_directory}/" if save_directory is not None else ""

    ingest_filename = f"{directory_prefix}{prefix}_ingest_info.avro"
    cell_filename = f"{directory_prefix}{prefix}_cell_info.avro"
    feature_filename = f"{directory_prefix}{prefix}_feature_info.avro"
    raw_counts_filename = f"{directory_prefix}{prefix}_raw_counts.avro"

    print(f"Loading input AnnData file '{input_file}'...")
    adata = optimized_read_andata(input_file)

    print("Processing ingest metadata...")
    dump_ingest_info(
        adata=adata,
        filename=ingest_filename,
        dataset_id=czi_dataset_id,
        dataset_version_id=czi_dataset_version_id,
        ingest_id=ingest_id,
        load_uns_data=load_uns_data,
        prefix=prefix,
        included_adata_uns_keys=included_adata_uns_keys,
        source_filename=input_file,
    )

    print("Processing cell/observation metadata...")
    dump_cell_info(adata=adata, filename=cell_filename, cas_cell_index_start=cas_cell_index_start, ingest_id=ingest_id)

    print("Processing feature/gene/variable metadata...")
    dump_feature_info(
        adata=adata,
        filename=feature_filename,
        cas_feature_index_start=cas_feature_index_start,
        ingest_id=ingest_id,
        original_feature_id_lookup=original_feature_id_lookup,
    )

    print("Processing core data...")

    total_cells = get_x_matrix(adata).shape[0]

    # close and attempt to free memory
    adata.file.close()
    del adata
    gc.collect()

    num_batches = math.ceil(total_cells / count_matrix_multiprocessing_batch_size)

    # ranges are start-inclusive and end-exclusive
    batches = [
        (
            x,
            x * count_matrix_multiprocessing_batch_size,
            min(total_cells, x * count_matrix_multiprocessing_batch_size + count_matrix_multiprocessing_batch_size),
        )
        for x in range(num_batches)
    ]

    start = current_milli_time()
    with multiprocessing.get_context("spawn").Pool() as pool:
        args = []
        for index, row_offset, end in batches:
            chunk_raw_counts_filename = raw_counts_filename.replace(".avro", f".{index:06}.csv")
            args.append(
                [input_file, row_offset, end, cas_cell_index_start, cas_feature_index_start, chunk_raw_counts_filename]
            )

        for result in pool.map(process_dump_core_matrix, args, chunksize=1):
            print(f"Got result: {result} for {args}", flush=True)

    print(f"    Processed {total_cells} cells... in {current_milli_time() - start} ms")

    print("Done.")


def process_dump_core_matrix(args):
    start = current_milli_time()

    (input_file, row_offset, end, cas_cell_index_start, cas_feature_index_start, raw_counts_filename) = args
    coord = optimized_read_raw_X(input_file, row_offset, end)
    print(f"    {raw_counts_filename} - Read anndata, subset matrix... in {current_milli_time() - start} ms")

    dump_core_matrix(
        coord.row, coord.col, coord.data, row_offset, cas_cell_index_start, cas_feature_index_start, raw_counts_filename
    )


# based on implementation of
# https://github.com/scverse/anndata/blob/6473f2034aa6e28ebc826ceeab15f413b8d294d8/anndata/_io/h5ad.py#L119
# optimized for reading of subset of raw.X
def optimized_read_andata(input_file):
    f = h5py.File(input_file, "r")

    attributes = ["obs", "var", "uns"]
    d = dict(filename=input_file, filemode="r")
    d.update({k: read_elem(f[k]) for k in attributes if k in f})

    d["raw"] = _read_raw(f, attrs={"var", "varm"})

    # Backwards compat to <0.7
    if isinstance(f["obs"], h5py.Dataset):
        _clean_uns(d)

    return AnnData(**d)


def optimized_read_raw_X(input_file, row_offset, end):
    adata = optimized_read_andata(input_file)
    coord = get_x_matrix(adata)[row_offset:end, :].tocoo()
    adata.file.close()
    del adata
    gc.collect()
    return coord
