#https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
#https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
#https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad
import gzip
import time
import sys
import argparse

def current_milli_time():
    return round(time.time() * 1000)

def dump_core_matrix(x, row_lookup, col_lookup):
    counter = 0
    start = current_milli_time()
    with gzip.open('cas_raw_counts.tsv.gz', 'wt') as f:
        f.write("cas_cell_index\tcas_gene_index\tdata\n")
        
        cx = x.tocoo()    
        for i,j,v in zip(cx.row, cx.col, cx.data):
            cas_cell_index = row_lookup[i]
            cas_gene_index = col_lookup[j]
            
            # Todo -- how can we ensure this is safe/right?
            v_int = int(v)
            f.write(f"{cas_cell_index}\t{cas_gene_index}\t{v_int}\n")
            counter = counter + 1
            if (counter % 1000000 == 0):
                end = current_milli_time()
                print(f"    Processed {counter} rows... in {end-start} ms")
                start = end

def process(input_file, cas_cell_index_start, cas_gene_index_start):
    print("Loading data...")
    adata = ad.read(input_file)
    
    # dump out cell info (obs) -- cell_id is index (26 columns)
    print("Processing observation metadata...")
    adata.obs['original_cell_id'] = adata.obs.index    
    adata.obs['cas_cell_index'] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    
    row_index_to_cas_cell_index = [None] * len(adata.obs)
    for row in range(0, len(adata.obs)):
        row_index_to_cas_cell_index[row] = adata.obs['cas_cell_index'].iloc[[row]][0]
    
    adata.obs[['cas_cell_index','original_cell_id','cell_type']].to_csv('cas_cell_info.tsv.gz', index=False, sep='\t', compression='gzip')
    
    # dump out gene info -- gene_id is index (12 columns)
    print("Processing variable metadata...")
    adata.var['original_gene_id'] = adata.var.index
    
    adata.var['cas_gene_index'] = np.arange(cas_gene_index_start, cas_gene_index_start + len(adata.var))
    col_index_to_cas_gene_index = [None] * len(adata.var)
    for col in range(0, len(adata.var)):
        col_index_to_cas_gene_index[col] = adata.var['cas_gene_index'].iloc[[row]][0]
    
    adata.var[['cas_gene_index','original_gene_id','feature_name']].to_csv('cas_gene_info.tsv.gz', index=False, sep='\t', compression='gzip')
    
    print("Processing core data...")
    
    ## recode the indexes to be the indexes of obs/var or should obs/var include these indicies?
    dump_core_matrix(adata.raw.X, row_index_to_cas_cell_index, col_index_to_cas_gene_index)
        
# TODO: naieve dump, compare
#    rows,cols = adata.X.nonzero()
#    for row,col in zip(rows,cols):
#        cas_cell_index = row_index_to_cas_cell_index[row]
#        cas_gene_index = col_index_to_cas_gene_index[col]
#        v = adata.X[row,col]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Convert AnnData Single Cell Expression Data into format for loading into BQ')

    parser.add_argument('--input',type=str, help='AnnData format input file', required=True)
    parser.add_argument('--cas_cell_index_start',type=int, help='starting number for cell index', required=False, default=0)
    parser.add_argument('--cas_gene_index_start',type=int, help='starting number for gene index', required=False, default=0)

    # Execute the parse_args() method
    args = parser.parse_args()

    process(args.input, args.cas_cell_index_start, args.cas_gene_index_start)



