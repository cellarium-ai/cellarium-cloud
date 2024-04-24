# Purpose

This is a python script that can be run from a single VM (provided the right credentials are loaded in .env)

The script will look at a specific google bucket folder, find all the h5ad files, and sequentially run avro creation and avro ingest to BigQuery for each file.

# Usage

Make a local working dirtory (will get populated with files):
```console
cd <repo_root>
mkdir ingest
cd ingest
```

Run the script (good idea to use a tmux session):
```console
python ../src/casp/workflows/local/ingest_from_bucket.py --dataset my_datastore --gcs_bucket_name my-cellarium-ingest-bucket --gcs_h5ad_dir my_h5ad_folder --gcs_avro_dir ingest_staging
```

- `my_datastore` is the (desired) name of the BigQuery dataset (need not exist, but it can exist with data)
- `my-cellarium-ingest-bucket` is the name of the google bucket that contains the h5ad files
- `my_h5ad_folder`: h5ad files are in `gs://my-cellarium-ingest-bucket/my_h5ad_folder/`
- `ingest_staging`: the folder `gs://my-cellarium-ingest-bucket/ingest_staging/` is created and deleted over and over for tmp files
