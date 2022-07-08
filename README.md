# Cell Annotation Service Pilot

This repository is intended to hold the functional prototype for a single cell annotation service

### Prerequisites / Installation

 - Python 3.7+
 
### Testing Data

For local testing, we are using a small data set from CellXGene.

[Horizontal Cells in Human Retina](https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7) - 7,348 cells

And a larger data set

[Tabula Sapiens Endothelial](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) -  31,691 cells

### Initialize BigQuery Data Set

This will create the dataset, and the core data model.  If any of these exist it will continue on

```
python src/initialize_dataset.py --project <google-project> --dataset <dataset-name>
```

### Convert AnnData to Load Format

For example, to load the retina data set

```
python src/anndata_to_bq.py --input data/horizontal-cells-in-human-retina.h5ad --cas_cell_index_start 1000 --cas_gene_index_start 5000
```

NOTE: if multiple data sets are being loaded into the same dataset, offsets for the cell and gene indexes must be provided via the `cas_cell_index_start` and `cas_gene_index_start` parameters respectively

### Ingest into BigQuery

Example using a specific projct (broad-dsp-spec-ops) and dataset (kc_cas_test_v1)

```
bq load -project_id broad-dsp-spec-ops -F tab --skip_leading_rows 1 kc_cas_test_v1.cas_cell_info cas_cell_info.tsv.gz
bq load -project_id broad-dsp-spec-ops -F tab --skip_leading_rows 1 kc_cas_test_v1.cas_gene_info cas_gene_info.tsv.gz
bq load -project_id broad-dsp-spec-ops -F tab --skip_leading_rows 1 kc_cas_test_v1.cas_raw_count_matrix cas_raw_counts.tsv.gz
```

### Extract Random Subset

This will randomly get a specified number of cells' data from your BigQuery dataset

```
python src/random_bq_to_anndata  --project <google-project> --dataset <dataset-name> --num_cells <num>
```

TODO: use data returned from BigQuery to output an AnnData file.
