# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

import anndata as ad
import argparse
from fastavro import writer, parse_schema
import hashlib
import numpy as np
import os
import time


def current_milli_time():
    return round(time.time() * 1000)


def write_avro(generator, parsed_schema, filename, progress_batch_size=10000):
    start = current_milli_time()
    with open(filename, 'a+b') as out:
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


def dump_core_matrix(x, row_lookup, col_lookup, filename):
    schema = {
        'doc': 'Raw data indexed by cell and feature id in the CAS BigQuery schema',
        'namespace': 'cas',
        'name': 'CellFeature',
        'type': 'record',
        'fields': [
            {'name': 'cas_cell_index', 'type': 'int'},
            {'name': 'cas_feature_index', 'type': 'int'},
            {'name': 'raw_counts', 'type': 'int'}
        ]
    }
    parsed_schema = parse_schema(schema)

    def raw_counts_generator():
        cx = x.tocoo(copy=False)

        for i, j, v in zip(cx.row, cx.col, cx.data):
            cas_cell_index = row_lookup[i]
            cas_feature_index = col_lookup[j]
            
            # Todo -- how can we ensure this is safe/right?
            v_int = int(v)
            yield {
                u'cas_cell_index': cas_cell_index.item(),
                u'cas_feature_index': cas_feature_index.item(),
                u'raw_counts': v_int
            }

    write_avro(raw_counts_generator, parsed_schema, filename, progress_batch_size=1000000)


def dump_cell_info(adata, filename, cas_cell_index_start, ingest_id):
    adata.obs['original_cell_id'] = adata.obs.index
    adata.obs['cas_cell_index'] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    row_index_to_cas_cell_index = [None] * len(adata.obs)
    for row in range(0, len(adata.obs)):
        row_index_to_cas_cell_index[row] = adata.obs['cas_cell_index'].iloc[[row]][0]

    schema = {
        'doc': 'CAS Cell Info',
        'namespace': 'cas',
        'name': 'CellInfo',
        'type': 'record',
        'fields': [
            {'name': 'cas_cell_index', 'type': 'int'},
            {'name': 'original_cell_id', 'type': 'string'},
            {'name': 'cell_type', 'type': 'string'},
            {'name': 'obs_metadata', 'type': 'string'},
            {'name': 'cas_ingest_id', 'type': 'string'},
        ]
    }

    parsed_schema = parse_schema(schema)
    added_columns = ['cas_cell_index', 'original_cell_id']
    obs_original = adata.obs.drop(columns=added_columns)
    cells = adata.obs[added_columns + ['cell_type']]

    def cell_generator():
        for i, r in cells.iterrows():
            yield {
                u'cas_cell_index': r['cas_cell_index'],
                u'original_cell_id': r['original_cell_id'],
                u'cell_type': r['cell_type'],
                u'obs_metadata': obs_original.loc[i].to_json(),
                u'cas_ingest_id': ingest_id,
            }

    write_avro(cell_generator, parsed_schema, filename)

    return row_index_to_cas_cell_index


def dump_feature_info(adata, filename, cas_feature_index_start, ingest_id):
    adata.var['original_feature_id'] = adata.var.index
    adata.var['cas_feature_index'] = np.arange(cas_feature_index_start, cas_feature_index_start + len(adata.var))

    col_index_to_cas_feature_index = [None] * len(adata.var)
    for col in range(0, len(adata.var)):
        col_index_to_cas_feature_index[col] = adata.var['cas_feature_index'].iloc[[col]][0]

    schema = {
        'doc': 'CAS Feature Info',
        'namespace': 'cas',
        'name': 'FeatureInfo',
        'type': 'record',
        'fields': [
            {'name': 'cas_feature_index', 'type': 'int'},
            {'name': 'original_feature_id', 'type': 'string'},
            {'name': 'feature_name', 'type': 'string'},
            {'name': 'var_metadata', 'type': 'string'},
            {'name': 'cas_ingest_id', 'type': 'string'},
        ]
    }

    parsed_schema = parse_schema(schema)
    added_columns = ['cas_feature_index', 'original_feature_id']
    var_original = adata.var.drop(columns=added_columns)
    features = adata.var[added_columns + ['feature_name']]

    def feature_generator():
        for i, row in features.iterrows():
            yield {
                u'cas_feature_index': row['cas_feature_index'],
                u'original_feature_id': row['original_feature_id'],
                u'feature_name': row['feature_name'],
                u'var_metadata': var_original.loc[i].to_json(),
                u'cas_ingest_id': ingest_id,
            }

    write_avro(feature_generator, parsed_schema, filename)
    return col_index_to_cas_feature_index


def confirm_output_files_do_not_exist(filenames):
    existing = list(filter(lambda f: os.path.exists(f), filenames))
    if len(existing) > 0:
        raise ValueError(
            f"Found existing output files, please rename or remove before running conversion: {', '.join(existing)}")


# https://stackoverflow.com/a/3431838
def md5(filename):
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def process(input_file, cas_cell_index_start, cas_feature_index_start, avro_prefix):
    avro_prefix = "cas" if not avro_prefix else avro_prefix
    print(f"Hashing input AnnData file '{input_file}'...")
    # This will look like 'cas-ingest-00d00dad'
    ingest_id = f'cas-ingest-{md5(input_file)[:8]}'

    file_types = ['cell_info', 'feature_info', 'raw_counts']
    filenames = [f'{avro_prefix}_{file_type}.avro' for file_type in file_types]
    confirm_output_files_do_not_exist(filenames)
    (cell_filename, feature_filename, raw_counts_filename) = filenames

    print("Loading data...")
    adata = ad.read(input_file)

    print("Processing cell/observation metadata...")
    row_index_to_cas_cell_index = dump_cell_info(adata, cell_filename, cas_cell_index_start, ingest_id)
    # dump out feature info -- feature_id is index (12 columns)
    print("Processing feature/gene/variable metadata...")
    col_index_to_cas_feature_index = dump_feature_info(adata, feature_filename, cas_feature_index_start, ingest_id)

    print("Processing core data...")
    # recode the indexes to be the indexes of obs/var or should obs/var include these indices?
    dump_core_matrix(adata.raw.X, row_index_to_cas_cell_index, col_index_to_cas_feature_index, raw_counts_filename)
        
# TODO: naive dump, compare
#    rows,cols = adata.X.nonzero()
#    for row,col in zip(rows,cols):
#        cas_cell_index = row_index_to_cas_cell_index[row]
#        cas_feature_index = col_index_to_cas_feature_index[col]
#        v = adata.X[row,col]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='Convert AnnData Single Cell Expression Data into format for loading into BQ')

    parser.add_argument('--input', type=str, help='AnnData format input file', required=True)
    parser.add_argument('--avro_prefix', type=str,
                        help='Prefix to use for output Avro files, e.g. <avro_prefix>_cell_info.avro etc',
                        required=False)
    parser.add_argument('--cas_cell_index_start', type=int, help='starting number for cell index', required=False,
                        default=0)
    parser.add_argument('--cas_feature_index_start', type=int, help='starting number for feature index', required=False,
                        default=0)
    parser.add_argument('--flush_batch_size', type=int, help='max size of Avro batches to flush', required=False)

    args = parser.parse_args()

    global flush_batch_size

    flush_batch_size = 10000 if not args.flush_batch_size else args.flush_batch_size

    process(args.input, args.cas_cell_index_start, args.cas_feature_index_start, args.avro_prefix)
