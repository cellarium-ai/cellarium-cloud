from kfp import dsl

from casp.workflows.kubeflow import machine_specs


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def create_avro_files(gcs_config_path: str):
    import os
    import tempfile

    import yaml
    from smart_open import open

    from casp.bq_scripts.anndata_to_avro import process as anndata_to_avro
    from casp.services import utils

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    gcs_bucket_name = config_data["gcs_bucket_name"]
    gcs_file_path = config_data["gcs_file_path"]
    uns_meta_keys = config_data["uns_meta_keys"]
    cas_cell_index_start = config_data["cas_cell_index_start"]
    cas_feature_index_start = config_data["cas_feature_index_start"]
    load_uns_data = config_data["load_uns_data"]
    original_feature_id_lookup = config_data["original_feature_id_lookup"]
    dataset_id = config_data["dataset_id"]
    dataset_version_id = config_data["dataset_version_id"]
    gcs_stage_dir = config_data["gcs_stage_dir"]

    temp_dir = tempfile.TemporaryDirectory()

    file_name = gcs_file_path.split("/")[-1]

    utils.download_file_from_bucket(
        bucket_name=gcs_bucket_name, source_blob_name=gcs_file_path, destination_file_name=file_name
    )
    _, project_id = utils.get_google_service_credentials()

    if uns_meta_keys is not None:
        uns_meta_keys_list = set([x.strip() for x in uns_meta_keys.split(",")])
    else:
        uns_meta_keys_list = None

    prefix = file_name.split(".")[0]

    anndata_to_avro(
        input_file=file_name,
        project=project_id,
        cas_cell_index_start=cas_cell_index_start,
        cas_feature_index_start=cas_feature_index_start,
        prefix=prefix,
        czi_dataset_id=dataset_id,
        czi_dataset_version_id=dataset_version_id,
        load_uns_data=load_uns_data,
        original_feature_id_lookup=original_feature_id_lookup,
        included_adata_uns_keys=uns_meta_keys_list,
        save_directory=temp_dir.name,
    )
    ingest_files = [x for x in os.listdir(temp_dir.name) if x.startswith(prefix) and not x.endswith(".h5ad")]

    for ingest_file_name in ingest_files:
        utils.upload_file_to_bucket(
            local_file_name=f"{temp_dir.name}/{ingest_file_name}",
            bucket=gcs_bucket_name,
            blob_name=f"{gcs_stage_dir}/{ingest_file_name}",
        )

    temp_dir.cleanup()


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def ingest_data(gcs_config_path: str):
    import math
    import pathlib
    import time

    import yaml
    from google.api_core.exceptions import Forbidden
    from google.cloud import bigquery
    from smart_open import open

    from casp.bq_scripts import create_bigquery_objects, ingest_data_to_bq
    from casp.services import utils

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    project_id = config_data["project_id"]
    gcs_bucket_name = config_data["gcs_bucket_name"]
    dataset = config_data["dataset"]
    gcs_stage_dir = config_data["gcs_stage_dir"]

    bq_client = bigquery.Client()
    create_bigquery_objects(client=bq_client, project=project_id, dataset=dataset)

    ingest_file_blobs = utils.list_blobs(bucket_name=gcs_bucket_name, prefix=gcs_stage_dir)
    blob_names = [x.name for x in ingest_file_blobs]
    ingest_avro_prefixes = []

    for blob_name in blob_names:
        path = pathlib.Path(blob_name)
        blob_directory = str(path.parent)

        # Getting prefixes (everything before `_cell_info`, `_ingest_info`, `_feature_info` `_raw_counts`)
        suffixes = ["_cell_info", "_ingest_info", "_feature_info", "_raw_counts"]
        prefix = path.name

        for suffix in suffixes:
            prefix = prefix.split(suffix)[0]

        if blob_directory != gcs_stage_dir:
            continue

        ingest_avro_prefixes.append(prefix)

    ingest_avro_prefixes = list(set(ingest_avro_prefixes))

    for avro_prefix in ingest_avro_prefixes:
        need_retry = True
        attempt_counter = 1

        while need_retry and attempt_counter <= 5:
            try:
                print(f"Ingesting files with prefix: {avro_prefix}")
                ingest_data_to_bq(
                    project=project_id,
                    dataset=dataset,
                    gcs_bucket_name=gcs_bucket_name,
                    avro_prefix=avro_prefix,
                    gcs_stage_dir=gcs_stage_dir,
                )
            except Forbidden as e:
                # It can happen when limit of number of table update operations is exceeded
                # Retrying the operation 5 times with an exponential backoff as suggested in the docs:
                # https://cloud.google.com/bigquery/quotas#standard_tables
                print("Was not able to ingest data", e)
                time_to_wait = math.exp(attempt_counter)
                time.sleep(time_to_wait)
                if attempt_counter <= 5:
                    print("Retrying another attempt...")
            except Exception as e:
                import pandas as pd

                bad_prefixes_name = f"gs://{gcs_bucket_name}/dev/ingest/bad_prefixes_missing.csv"
                try:
                    df = pd.read_csv(bad_prefixes_name)
                except FileNotFoundError:
                    df = pd.DataFrame(columns=("name_prefix", "exception_text"))

                new_row = pd.DataFrame({"name_prefix": [avro_prefix], "exception_text": [str(e)]})
                df = pd.concat([df, new_row], ignore_index=True)
                df.to_csv(bad_prefixes_name, index=False)
                need_retry = False

            else:
                need_retry = False

            attempt_counter += 1


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def precalculate_fields(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.bq_scripts import precalculate_fields

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    fields_str = config_data["fields"]
    dataset = config_data["dataset"]
    project_id = config_data["project_id"]

    fields_list = fields_str.split(",")
    precalculate_fields(dataset=dataset, fields=fields_list, project=project_id)


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def prepare_extract(gcs_config_path: str):
    import json
    import tempfile

    import yaml
    from smart_open import open

    from casp.bq_scripts import prepare_all_cell_types, prepare_extract, prepare_measured_genes_info
    from casp.services import utils

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    filters_json_path = config_data["filters_json_path"]
    obs_columns_to_include = config_data["obs_columns_to_include"]
    project_id = config_data["project_id"]
    dataset = config_data["dataset"]
    extract_table_prefix = config_data["extract_table_prefix"]
    fq_allowed_original_feature_ids = config_data["fq_allowed_original_feature_ids"]
    extract_bin_size = config_data["extract_bin_size"]
    bucket_name = config_data["bucket_name"]
    extract_bucket_path = config_data["extract_bucket_path"]

    obs_columns_to_include_split = None if obs_columns_to_include is None else obs_columns_to_include.split(",")

    filters = None

    with open(filters_json_path) as f:
        filters = json.loads(f.read())

    prepare_extract(
        project=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        filters=filters,
        obs_columns_to_include=obs_columns_to_include_split,
    )
    measured_genes_info_df = prepare_measured_genes_info(
        project=project_id,
        dataset=dataset,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
    )
    all_cell_types_df = prepare_all_cell_types(project=project_id, dataset=dataset)

    with tempfile.TemporaryDirectory() as temp_dir:
        measured_genes_file_name = f"measured_genes_info.csv"
        measured_genes_file_path = f"{temp_dir}/{measured_genes_file_name}"
        measured_genes_info_df.to_csv(measured_genes_file_path)

        all_cell_types_file_name = f"all_cell_types.csv"
        all_cell_types_file_path = f"{temp_dir}/{all_cell_types_file_name}"
        all_cell_types_df.to_csv(all_cell_types_file_path, index=False)

        utils.upload_file_to_bucket(
            local_file_name=measured_genes_file_path,
            bucket=bucket_name,
            blob_name=f"{extract_bucket_path}/shared_meta/{measured_genes_file_name}",
        )
        utils.upload_file_to_bucket(
            local_file_name=all_cell_types_file_path,
            bucket=bucket_name,
            blob_name=f"{extract_bucket_path}/shared_meta/{all_cell_types_file_name}",
        )


@dsl.component(base_image=machine_specs.DOCKER_IMAGE_NAME_CPU)
def extract(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.bq_scripts import extract_bins_in_parallel_workers

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    obs_columns_to_include = config_data["obs_columns_to_include"]
    project_id = config_data["project_id"]
    dataset = config_data["dataset"]
    extract_table_prefix = config_data["extract_table_prefix"]
    output_bucket_name = config_data["output_bucket_name"]
    extract_bucket_path = config_data["extract_bucket_path"]
    start_bin = config_data["start_bin"]
    end_bin = config_data["end_bin"]
    obs_columns_to_include_list = obs_columns_to_include.split(",")

    extract_bins_in_parallel_workers(
        project_id=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=start_bin,
        end_bin=end_bin,
        output_bucket_name=output_bucket_name,
        extract_bucket_path=extract_bucket_path,
        obs_columns_to_include=obs_columns_to_include_list,
    )
