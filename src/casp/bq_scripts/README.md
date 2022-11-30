## BQ Data Operation Scripts

### Convert AnnData to Avro Load Format

First convert from AnnData to the Avro format that will be used to load BigQuery. The `anndata_to_avro.py` script needs
to know the appropriate values to use for cell and feature indexes. The script can either be pointed to an existing
dataset / project to figure those out for itself, or these values can be be specified as command line arguments.

To figure out start values from a dataset / project:

```
python src/casp/bq_scipts/anndata_to_avro.py --input data/horizontal-cells-in-human-retina.h5ad --avro_prefix retina --project <google-project> --dataset <dataset-name>
```

As specified command line arguments:

```
python src/casp/bq_scipts/anndata_to_avro.py --input data/horizontal-cells-in-human-retina.h5ad --avro_prefix retina --cas_cell_index_start 1000 --cas_feature_index_start 5000
```

### Initialize dataset and ingest Avro files into BigQuery

Continuing this example:

```
python src/casp/bq_scipts/load_dataset.py --project <google-project> --dataset <dataset-name> --avro_prefix retina --gcs_path_prefix <gcs-prefix-to-avro-files>
```

All of these parameters are required.

* `project` specifies the Google project id that will contain the BigQuery dataset.
* `dataset` specifies the dataset name that will contain the BigQuery tables.
* `avro_prefix` specifies the prefix used to name the output Avro from the AnnData conversion in the previous step.
* `gcs_path_prefix` specifies the GCS prefix to which the Avro files should be staged for BigQuery ingestion.

This step:

* Creates the dataset in this specified project if required
* Creates the CASP tables in the specified dataset if required
* Stages the specified Avro files to the specified GCS path
* Loads the specified tables from the staged Avro files
* Cleans up the staged Avro files from GCS


## Data Extract

Extracting data into minibatchs of N cells is comprised of two steps:

* Prepare - randomizes, preprocesses and shards the entire dataset
* Extract - extracts one or more shards, can be done in parallel

### Prepare

To prepare a dataset, run:

```
python src/casp/bq_scipts/prepare_for_training.py --project <project> --dataset <dataset> --extract_table_prefix <table_prefix> --extract_bin_size <bin_size>
```

where

* `project` specifies the Google project id that will contain the BigQuery dataset.
* `dataset` specifies the dataset name that will contain the BigQuery tables.
* `table_prefix` specifies the prefix used to name the temporary tables
* `bin_size` approximate number of cells for each shard/bin (e.g. 10000)

### Extract

```
python src/casp/bq_scipts/extract_minibatch_to_anndata.py --project <project> --dataset <dataset> --extract_table_prefix <table_prefix> --start_bin <start_bin> --end_bin <end_bin> --output <output_file_name>
```

where

* `project` specifies the Google project id that will contain the BigQuery dataset.
* `dataset` specifies the dataset name that will contain the BigQuery tables.
* `table_prefix` specifies the prefix used to name the temporary tables
* `start_bin` 0-based, inclusive lower bound bin number for this extract
* `end_bin` 0-based, inclusive upper bound bin number for this extract (can be the same as start for single bin)
* `output_file_name` output file (e.g.`demo_4m_v2.0000.h5ad`)

### Extract Random Subset without Prepare (deprecated)

This will randomly get a specified number of cells' data from your BigQuery dataset:

```
python src/casp/bq_scipts/random_bq_to_anndata.py --project <google-project> --dataset <dataset-name> --num_cells <num> --output_file_prefix <output prefix>
```

One AnnData file will be written for each ingest (input AnnData file) in which the randomly selected cells were loaded.
e.g. if 100 cells are randomly selected that were ingested from three separate AnnData files, three output files will be
written. Files will have names like `<output prefix>-<ingest id>.h5ad`. Each file will contain only the features
associated with that particular ingest, and only the sub-subset of cells from that ingest within the random subset.
Features and cells will be assigned `var` and `obs` metadata respectively, and the overall AnnData file will have its
original `uns` metadata.