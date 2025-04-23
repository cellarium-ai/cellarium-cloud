# Purpose

This is a python script that can be run from a single VM (provided the right credentials are loaded in .env)

The script will look at a specific google bucket folder, find all the h5ad files, and sequentially run avro creation and avro ingest to BigQuery for each file.

# Usage

Ensure a `.env` file exists at `cellarium-cloud/src/casp/services/.env` and is populated appropriately with the result of downloading a service account key from a Google service account (created for this purpose) on the relevant Google project with `BigQuery Admin` permissions. Contents of that file should look like this

```
GOOGLE_SERVICE_ACCOUNT_CREDENTIALS={
    "type": "service_account",
    "project_id": "my-google-project",
    "private_key_id": "somenumbersandletters",
    "private_key": "-----BEGIN PRIVATE KEY-----\nbunchofnumbersandletters\n-----END PRIVATE KEY-----\n",
    "client_email": "named-service-account@my-google-project.iam.gserviceaccount.com",
    "client_id": "somenumber",
    "auth_uri": "https://accounts.google.com/o/oauth2/auth",
    "token_uri": "https://oauth2.googleapis.com/token",
    "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
    "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/named-service-account%40my-google-project.iam.gserviceaccount.com",
    "universe_domain": "googleapis.com"
}
```

Make a local working dirtory (will get populated with files):
```console
cd <repo_root>
mkdir ingest
cd ingest
```

Make the code visible
```console
export PYTHONPATH="${PYTHONPATH}:/path/to/cellarium-cloud/src"
```

Run the script (good idea to use a tmux session):
```console
python ../src/casp/workflows/local/ingest_from_bucket.py --dataset my_datastore --gcs_bucket_name my-cellarium-ingest-bucket --gcs_h5ad_dir my_h5ad_folder --gcs_avro_dir ingest_staging
```

- `my_datastore` is the (desired) name of the BigQuery dataset (need not exist, but it can exist with data)
- `my-cellarium-ingest-bucket` is the name of the google bucket that contains the h5ad files
- `my_h5ad_folder`: h5ad files are in `gs://my-cellarium-ingest-bucket/my_h5ad_folder/`
- `ingest_staging`: the folder `gs://my-cellarium-ingest-bucket/ingest_staging/` is created and deleted over and over for tmp files
