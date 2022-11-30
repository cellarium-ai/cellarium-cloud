# Cell Annotation Service Pilot

This repository is intended to hold the functional prototype for a single cell annotation service

### Prerequisites / Installation

 - Python 3.7+

### Developer Setup

Some #protips:

* Ingest will often time out when running from a local machine, which will leave the BigQuery dataset and tables in a weird state. It's probably a good idea to spin up a GCE VM if ingesting non-toy datasets.
* Some of the larger AnnData files from cellxgene can require a lot of memory for Avro conversion (e.g. > 64 GiB for the full trisomy 18 dataset); another reason to spin up a GCE VM for non-toy ingests.

To create a virtual python environment:

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




## Reference Data Preparation

To extract data across multiple datasets, we have adopted an "allow list" of ensemble genes based on 10X Genomics reference data.

```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
tail -n+2 refdata-gex-GRCh38-2020-A/star/geneInfo.tab > /tmp/gene.tsv
bq load dsp-cell-annotation-service:cas_reference_data.refdata-gex-GRCh38-2020-A /tmp/gene.tsv feature_name:STRING
```
