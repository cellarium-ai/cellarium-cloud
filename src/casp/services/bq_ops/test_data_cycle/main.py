import os
import typing as t

import anndata
from google.cloud import bigquery

from casp.services import utils
from casp.services.bq_ops.anndata_to_ingest_files.main import main as anndata_to_ingest_files
from casp.services.bq_ops.extract.main import main as extract
from casp.services.bq_ops.ingest_files_to_bq.main import main as ingest_files_to_bq
from casp.services.bq_ops.prepare_extract.main import main as prepare_extract
from casp.services.bq_ops.test_data_cycle import constants, data_controller

successes = []
failures = []


def get_source_file(gcs_bucket_name: str, gcs_input_file_paths: t.List[str]) -> "anndata.AnnData":
    input_file_names_local = []

    for i, input_file_path in enumerate(gcs_input_file_paths):
        input_file_name = f"input_file_{i}.h5ad"
        input_file_names_local.append(input_file_name)
        utils.download_file_from_bucket(
            bucket_name=gcs_bucket_name, source_blob_name=input_file_path, destination_file_name=input_file_name
        )

    adatas = [anndata.read_h5ad(x) for x in input_file_names_local]

    print("Cleaning up input files from disk")
    for f_name in input_file_names_local:
        os.remove(f_name)

    return anndata.concat(adatas=adatas, join="outer")


def get_anndata_or_raw(adata):
    if adata.raw is not None:
        return adata.raw
    else:
        return adata


def clean_up_cloud_from_test_case(bucket_name: str, dataset: str, extract_table_prefix: str):
    test_extract_data_dir = f"{extract_table_prefix}__data"
    dataset_info_dir = f"{dataset}_{extract_table_prefix}_info"
    utils.delete_folder_from_bucket(bucket_name=bucket_name, folder_name=test_extract_data_dir)
    utils.delete_folder_from_bucket(bucket_name=bucket_name, folder_name=dataset_info_dir)


def clean_up_cloud_from_tests_all():
    utils.delete_folder_from_bucket(bucket_name=constants.GCS_BUCKET_NAME, folder_name=constants.GCS_STAGE_DIR)
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    client.delete_dataset(constants.DATASET_NAME, delete_contents=True, not_found_ok=True)


def compare_extract_with_input_files(
    gcs_bucket_name: str,
    gcs_input_file_paths: t.List[str],
    gcs_extract_directory: str,
    dataset_name: str,
    number_of_extract_bins: int,
    test_id: str,
):
    extract_file_names = [f"extract_{x}.h5ad" for x in range(number_of_extract_bins)]
    adata_source = get_source_file(gcs_bucket_name=gcs_bucket_name, gcs_input_file_paths=gcs_input_file_paths)
    extracts_correspond_to_source = []

    for extract_file_name in extract_file_names:
        extract_chunk_name = f"{gcs_extract_directory}/{extract_file_name}"
        utils.download_file_from_bucket(
            bucket_name=gcs_bucket_name, source_blob_name=extract_chunk_name, destination_file_name=extract_file_name
        )
        adata_extract_chunk = anndata.read_h5ad(extract_file_name)
        # We take intersect features because some of the features that are present in extract files
        # weren't present in source anndata file
        intersect_features = list(
            set(adata_extract_chunk.var.index.values).intersection(set(adata_source.var.index.values))
        )

        cas_cell_ids = adata_extract_chunk.obs.index.values.tolist()

        original_cell_ids_map = dict(
            data_controller.get_original_cell_ids_by(dataset_name=dataset_name, cas_cell_ids=cas_cell_ids)
        )
        original_cell_ids = [original_cell_ids_map[x] for x in cas_cell_ids]

        # count_matrix_chunk = adata_extract_chunk.raw[:, intersect_features].X
        count_matrix_chunk = get_anndata_or_raw(adata_extract_chunk)[:, intersect_features].X
        count_matrix_source_file = get_anndata_or_raw(adata_source)[original_cell_ids, intersect_features].X

        values_are_same = (count_matrix_chunk != count_matrix_source_file).sum() == 0
        extracts_correspond_to_source.append(values_are_same)
        os.remove(extract_file_name)

    if all(extracts_correspond_to_source):
        successes.append(test_id)
        print("Success!")
    else:
        failures.append(test_id)
        print(f"Fail: {sum(extracts_correspond_to_source)} chunks are wrong")


def set_up_test_bq_database():
    ingest_index = 0
    for input_file_path in constants.GCS_INPUT_FILE_PATHS_ALL:
        print(f"Creating ingest files for {input_file_path}...")
        anndata_to_ingest_files(
            gcs_bucket_name=constants.GCS_BUCKET_NAME,
            gcs_file_path=input_file_path,
            cas_cell_index_start=ingest_index,
            cas_feature_index_start=ingest_index,
            original_feature_id_lookup="index",
            gcs_stage_dir=constants.GCS_STAGE_DIR,
        )
        ingest_index += 2000000
    print("Ingesting files to BigQuery dataset...")
    ingest_files_to_bq(
        dataset=constants.DATASET_NAME,
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_stage_dir=constants.GCS_STAGE_DIR,
        delete_ingest_files=False,
    )


def test_extract_filtered_by_homo_sapiens():
    test_extract_data_dir = f"{constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX}__data"
    print("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Homo sapiens",
    )
    print("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=8,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
    )
    print("Testing results...")
    compare_extract_with_input_files(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_directory=test_extract_data_dir,
        dataset_name=constants.DATASET_NAME,
        number_of_extract_bins=9,
        test_id="test_extract_filtered_by_homo_sapiens",
    )
    print("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_mus_mus():
    test_extract_data_dir = f"{constants.MUS_MUS_EXTRACT_TABLE_PREFIX}__data"
    print("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.MUS_MUS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Mus musculus",
    )
    print("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=2,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
    )
    print("Testing results...")
    compare_extract_with_input_files(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_MUS_MUS,
        gcs_extract_directory=test_extract_data_dir,
        dataset_name=constants.DATASET_NAME,
        number_of_extract_bins=3,
        test_id="test_extract_filtered_by_mus_mus",
    )
    print("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_homo_sapiens_small_chunk_size():
    test_extract_data_dir = f"{constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX}__data"
    print("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=5000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Homo sapiens",
    )
    print("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=17,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
    )
    print("Testing results...")
    compare_extract_with_input_files(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_directory=test_extract_data_dir,
        dataset_name=constants.DATASET_NAME,
        number_of_extract_bins=18,
        test_id="test_extract_filtered_by_homo_sapiens_small_chunk_size",
    )
    print("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_datasets():
    test_extract_data_dir = f"{constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX}__data"
    print("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_datasets="tg66rtgh7y-test-source-file-1_ingest_info.avro,yt55tgy4o-test-source-file-3_ingest_info.avro",
    )
    print("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=5,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
    )
    print("Testing results...")
    compare_extract_with_input_files(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_directory=test_extract_data_dir,
        dataset_name=constants.DATASET_NAME,
        number_of_extract_bins=5,
        test_id="test_extract_filtered_by_datasets",
    )
    print("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX,
    )


def main():
    """
    Test whether CAS infrastructure produces relevant anndata file chunks after adding them to the database, harmonizing
    and shuffling.
    """
    # Set up BigQuery Database
    set_up_test_bq_database()
    # Process test cases
    test_extract_filtered_by_homo_sapiens()
    test_extract_filtered_by_mus_mus()
    test_extract_filtered_by_homo_sapiens_small_chunk_size()
    test_extract_filtered_by_datasets()

    success_messages = "\n".join([f"\t{x}" for x in successes])
    failure_messages = "\n".join([f"\t{x}" for x in failures])

    print("======================================")

    print(f"{len(successes)} tests passed: \n {success_messages}")
    if len(failures) > 0:
        print(f"{len(failures)} tests failed: \n {failure_messages}")

    success_rate = (len(successes) / (len(successes) + len(failures))) * 100
    print(f"======= {success_rate:.0f}% succeeded! =======")

    print("Cleaning Up Infrastructure...")
    clean_up_cloud_from_tests_all()

    if len(failures) > 0:
        # Raise error so that if the test performed from Cromwell pipeline
        # it would mark it as failure in cromwell list
        raise ValueError("Not all the tests succeeded")

    print("Done.")


if __name__ == "__main__":
    main()
