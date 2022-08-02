# Cell Annotation Service Pilot

This repository is intended to hold the functional prototype for a single cell annotation service

### Prerequisites / Installation

 - Python 3.7+

### Developer Setup

Create a virtual python environment:

```shell
    python3 -mvenv casp-venv
    source casp-venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt
    pip install -r dev-requirements.txt
    pip install -e .
    pip install tox
```

To run unit tests:

```shell
    tox -e unit
```

To lint:

```shell
    tox -e lint
```

To automatically fix formatting issues:
```shell
    tox -e format
```
 
### Testing Data
    
For local testing, we are using a small data set from CellXGene.

[Horizontal Cells in Human Retina](https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7) - 7,348 cells

And a larger data set

[Tabula Sapiens Endothelial](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) -  31,691 cells


### Convert AnnData to Avro Load Format

First convert from AnnData to the Avro format that will be used to load BigQuery. The `anndata_to_avro.py` script needs
to know the appropriate values to use for cell and feature indexes. The script can either be pointed to an existing 
dataset / project to figure those out for itself, or these values can be be specified as command line arguments.

To figure out start values from a dataset / project:

```
python src/casp/anndata_to_avro.py --input data/horizontal-cells-in-human-retina.h5ad --avro_prefix retina --project <google-project> --dataset <dataset-name>
```

As specified command line arguments:

```
python src/casp/anndata_to_avro.py --input data/horizontal-cells-in-human-retina.h5ad --avro_prefix retina --cas_cell_index_start 1000 --cas_feature_index_start 5000
```

### Initialize dataset and ingest Avro files into BigQuery

Continuing this example:

```
python src/casp/load_dataset.py --project <google-project> --dataset <dataset-name> --avro_prefix retina --gcs_prefix <gcs-prefix-to-avro-files>
```

All of these parameters are required.

* `project` specifies the Google project id that will contain the BigQuery dataset.
* `dataset` specifies the dataset name that will contain the BigQuery tables.
* `avro_prefix` specifies the prefix used to name the output Avro from the AnnData conversion in the previous step.
* `gcs_prefix` specifies the GCS prefix to which the Avro files should be staged for BigQuery ingestion.

This step:

* Creates the dataset in this specified project if required
* Creates the CASP tables in the specified dataset if required
* Stages the specified Avro files to the specified GCS path
* Loads the specified tables from the staged Avro files
* Cleans up the staged Avro files from GCS


### Extract Random Subset

This will randomly get a specified number of cells' data from your BigQuery dataset:

```
python src/casp/random_bq_to_anndata.py --project <google-project> --dataset <dataset-name> --num_cells <num> --output_file_prefix <output prefix>
```

One AnnData file will be written for each ingest (input AnnData file) in which the randomly selected cells were loaded.
e.g. if 100 cells are randomly selected that were ingested from three separate AnnData files, three output files will be
written. Files will have names like `<output prefix>-<ingest id>.h5ad`. Each file will contain only the features
associated with that particular ingest, and only the sub-subset of cells from that ingest within the random subset.
Features and cells will be assigned `var` and `obs` metadata respectively, and the overall AnnData file will have its
original `uns` metadata.
