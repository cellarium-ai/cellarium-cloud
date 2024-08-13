from kfp import dsl


# @dsl.component()
def create_avro_files(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts.bq_ops.anndata_to_avro import create_ingest_files

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

    if uns_meta_keys is not None:
        uns_meta_keys_list = set([x.strip() for x in uns_meta_keys.split(",")])
    else:
        uns_meta_keys_list = None

    create_ingest_files(
        gcs_bucket_name=gcs_bucket_name,
        gcs_file_path=gcs_file_path,
        uns_meta_keys=uns_meta_keys_list,
        cas_cell_index_start=cas_cell_index_start,
        cas_feature_index_start=cas_feature_index_start,
        load_uns_data=load_uns_data,
        original_feature_id_lookup=original_feature_id_lookup,
        dataset_id=dataset_id,
        dataset_version_id=dataset_version_id,
        gcs_stage_dir=gcs_stage_dir,
    )


# @dsl.component()
def ingest_data(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts.bq_ops import ingest_data_to_bq

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    project_id = config_data["project_id"]
    gcs_bucket_name = config_data["gcs_bucket_name"]
    dataset = config_data["dataset"]
    gcs_stage_dir = config_data["gcs_stage_dir"]
    gcs_error_file_path = config_data["gcs_error_file_path"]

    ingest_data_to_bq(
        project_id=project_id,
        gcs_bucket_name=gcs_bucket_name,
        dataset=dataset,
        gcs_stage_dir=gcs_stage_dir,
        gcs_error_file_path=gcs_error_file_path,
    )


# @dsl.component()
def precalculate_fields(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts.bq_ops import precalculate_fields

    with open(gcs_config_path, mode="r") as file:
        config_data = yaml.safe_load(file)

    fields_str = config_data["fields"]
    dataset = config_data["dataset"]
    project_id = config_data["project_id"]

    fields_list = fields_str.split(",")
    precalculate_fields(dataset=dataset, fields=fields_list, project=project_id)


# @dsl.component()
def prepare_extract(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts.bq_ops import prepare_extract

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

    prepare_extract(
        project_id=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        fq_allowed_original_feature_ids=fq_allowed_original_feature_ids,
        extract_bin_size=extract_bin_size,
        filters_json_path=filters_json_path,
        obs_columns_to_include=obs_columns_to_include,
        bucket_name=bucket_name,
        extract_bucket_path=extract_bucket_path,
    )


# @dsl.component()
def extract(gcs_config_path: str):
    import yaml
    from smart_open import open

    from casp.scripts.bq_ops import extract_bins_in_parallel_workers

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

    extract_bins_in_parallel_workers(
        project_id=project_id,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        start_bin=start_bin,
        end_bin=end_bin,
        output_bucket_name=output_bucket_name,
        extract_bucket_path=extract_bucket_path,
        obs_columns_to_include=obs_columns_to_include,
    )
