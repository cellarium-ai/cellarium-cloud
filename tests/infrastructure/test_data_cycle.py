"""
Test utilities for data extraction processes.

This module contains test utilities for validating data extraction processes, including various
configurations like filtering by organism or dataset. Tests emphasize the correct preparation,
extraction, verification, and cleanup of resources. Input files are specifically designed for these test cases with
varying numbers of features.

Each test function:
1. Prepares the necessary extract tables.
2. Extracts the relevant data.
3. Verifies the correctness of the extracted chunks.
4. Cleans up resources after testing.
"""

import logging
import tempfile
import typing as t

import anndata
import pytest
from google.cloud import bigquery, exceptions

from casp.services import utils
from casp.services.bq_ops.anndata_to_ingest_files.main import main as anndata_to_ingest_files
from casp.services.bq_ops.extract.main import main as extract
from casp.services.bq_ops.ingest_files_to_bq.main import main as ingest_files_to_bq
from casp.services.bq_ops.precalculate_fields.main import main as precalculate_fields
from casp.services.bq_ops.prepare_extract.main import main as prepare_extract

from . import constants, data_controller

logger = logging.getLogger(__name__)


def concatenate_source_files(gcs_bucket_name: str, gcs_input_file_paths: t.List[str]) -> "anndata.AnnData":
    adatas = []
    for i, input_file_path in enumerate(gcs_input_file_paths):
        with tempfile.TemporaryFile() as f:
            utils.write_to_file_from_bucket(bucket_name=gcs_bucket_name, source_blob_name=input_file_path, file=f)
            adatas.append(anndata.read_h5ad(f))

    return anndata.concat(adatas=adatas, join="outer")


def get_anndata_or_raw(adata: "anndata.AnnData") -> "anndata.AnnData":
    """Retrieve the :class:`anndata.AnnData` object or its ``raw`` attribute, if it exists."""
    if adata.raw is not None:
        return adata.raw
    else:
        return adata


def clean_up_cloud_from_test_case(bucket_name: str, extract_table_prefix: str):
    """Clean up GCS from files that were produced by a particular test case."""
    test_extract_data_dir = f"{extract_table_prefix}__data"
    dataset_info_dir = f"{constants.DATASET_NAME}_{extract_table_prefix}_info"
    utils.delete_folder_from_bucket(bucket_name=bucket_name, folder_name=test_extract_data_dir)
    utils.delete_folder_from_bucket(bucket_name=bucket_name, folder_name=dataset_info_dir)


def clean_up_bq_database_and_gcs():
    """Clean up GCS and BigQuery from sources that were produced by all the tests."""
    utils.delete_folder_from_bucket(bucket_name=constants.GCS_BUCKET_NAME, folder_name=constants.GCS_STAGE_DIR)
    credentials, project = utils.get_google_service_credentials()
    client = bigquery.Client(credentials=credentials, project=project)
    client.delete_dataset(constants.DATASET_NAME, delete_contents=True, not_found_ok=True)


def set_up_bq_database():
    """Set up test BigQuery database"""
    ingest_index = 0
    for input_file_path in constants.GCS_INPUT_FILE_PATHS_ALL:
        logger.info(f"Creating ingest files for {input_file_path}...")
        anndata_to_ingest_files(
            gcs_bucket_name=constants.GCS_BUCKET_NAME,
            gcs_file_path=input_file_path,
            cas_cell_index_start=ingest_index,
            cas_feature_index_start=ingest_index,
            original_feature_id_lookup="index",
            gcs_stage_dir=constants.GCS_STAGE_DIR,
        )
        # Large number to be sure the start index of the next file will not repeat the index from the previous one
        ingest_index += 2000000
    logger.info("Ingesting files to BigQuery dataset...")
    ingest_files_to_bq(
        dataset=constants.DATASET_NAME,
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_stage_dir=constants.GCS_STAGE_DIR,
        delete_ingest_files=False,
    )
    logger.info("Precalculating fields...")
    precalculate_fields(dataset=constants.DATASET_NAME, fields_str=constants.PRECALCULATE_FIELDS)


@pytest.fixture(scope="module", autouse=True)
def manage_test_bq_database():
    """
    Pytest fixture to set up the BigQuery test database.

    This fixture is automatically invoked for tests within the module scope. It performs the following operations:
    1. Logs the initiation of BigQuery setup for test cases.
    2. Iterates through GCS input file paths, creating and logging ingest files for each path.
    3. Ingests the files to the specified BigQuery dataset.
    4. Precalculates specific fields for the dataset.

    After tests within the module scope are executed, it cleans up the associated cloud infrastructure files.

    Note:
        This fixture has a module scope and is autouse, meaning it is automatically invoked for every test in the
        module without needing to explicitly call or parameterize it in the test function.

    Yields:
        None: Just serves as a pause point for the test execution until the teardown/cleanup actions need to be
            executed.

    Raises:
        Any exceptions raised by the called functions like `anndata_to_ingest_files`, `ingest_files_to_bq`, etc.
    """
    logger.info("Set up BigQuery for the test cases")
    set_up_bq_database()
    yield
    logger.info("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_bq_database_and_gcs()


def verify_chunk_count_matrix(
    adata_extract_chunk: "anndata.AnnData", adata_source: "anndata.AnnData", dataset_name: str
) -> None:
    """
    Compare the extracted chunk raw count matrix with the source file. Use intersecting features because some of the
    features present in the extracted files aren't present in the source anndata file.

    :param adata_extract_chunk: Extracted chunk for verification.
    :param adata_source: The source file to compare against.
    :param dataset_name: Name of the dataset, utilized to retrieve original cell IDs.
    :raises AssertionError:  if chunk's raw counts correspond to source counts
    """

    intersect_features = list(
        set(adata_extract_chunk.var.index.values).intersection(set(adata_source.var.index.values))
    )

    cas_cell_ids = adata_extract_chunk.obs.index.values.tolist()

    original_cell_ids_map = dict(
        data_controller.get_original_cell_ids_by(dataset_name=dataset_name, cas_cell_ids=cas_cell_ids)
    )
    original_cell_ids = [original_cell_ids_map[x] for x in cas_cell_ids]

    count_matrix_chunk = get_anndata_or_raw(adata_extract_chunk)[:, intersect_features].X
    count_matrix_source_file = get_anndata_or_raw(adata_source)[original_cell_ids, intersect_features].X

    assert (
        count_matrix_chunk != count_matrix_source_file
    ).sum() == 0, "Count matrix does not correspond to values from source files"


def verify_chunk_total_mrna_umis(adata_extract_chunk: "anndata.AnnData", dataset_name: str) -> None:
    """
    Verify the ``total_mrna_umis`` of an extracted chunk against counts from BigQuery.

    This function retrieves mRNA counts from BigQuery based on the provided dataset name and CAS cell indexes.
    These indexes are taken from :attr:`adata_extract_chunk`. The function then compares these raw counts
    to the ``obs["total_mrna_umis"]`` values in the provided extracted chunk.

    :param adata_extract_chunk: Extracted chunk of data for verification.
    :param dataset_name: Name of the dataset, used to retrieve and compare original cell IDs from BigQuery.
    :raises AssertionError: If the chunk's `total_mrna_umis` doesn't match the calculated mRNA counts from BigQuery
    """
    cas_cell_ids = adata_extract_chunk.obs.index.values.tolist()
    calculated_mrna_counts = data_controller.get_raw_mrna_counts_by(
        dataset_name=dataset_name, cas_cell_ids=cas_cell_ids
    )

    total_mrna_umis = adata_extract_chunk.obs.loc[
        [x[0] for x in calculated_mrna_counts], "total_mrna_umis"
    ].values.tolist()
    calculated_mrna_counts = [x[1] for x in calculated_mrna_counts]
    assert (
        total_mrna_umis == calculated_mrna_counts
    ), "Values in obs column `total_mrna_umis` do not correspond to total mRNA counts from BigQuery"


def verify_extracted_chunks(
    gcs_bucket_name: str,
    gcs_input_file_paths: t.List[str],
    gcs_extract_file_paths: t.List[str],
    dataset_name: str,
):
    """
    Verifies extracted chunks of data from the Google Cloud Storage (GCS) against the source data.

    :param gcs_bucket_name: Name of the GCS bucket.
    :param gcs_input_file_paths: List of input file paths within the GCS bucket to verify against.
    :param gcs_extract_file_paths: List of extract file paths within the GCS bucket to verify with.
    :param dataset_name: Name of the BigQuery dataset.
    :raises AssertionError: If not all the chunks are verified
    """
    adata_source = concatenate_source_files(gcs_bucket_name=gcs_bucket_name, gcs_input_file_paths=gcs_input_file_paths)

    for extract_chunk_name in gcs_extract_file_paths:
        with tempfile.TemporaryFile() as f:
            try:
                utils.write_to_file_from_bucket(
                    bucket_name=gcs_bucket_name, source_blob_name=extract_chunk_name, file=f
                )
            except exceptions.NotFound:
                raise AssertionError(
                    f"Couldn't find an extract file `{extract_chunk_name}` in the GCS bucket. Please, make sure you "
                    f"called extract script before running this verification in your test case. If you did run run it, "
                    f"the error is probably in there."
                )

            adata_extract_chunk = anndata.read_h5ad(f)

            verify_chunk_count_matrix(
                adata_extract_chunk=adata_extract_chunk, adata_source=adata_source, dataset_name=dataset_name
            )
            verify_chunk_total_mrna_umis(adata_extract_chunk=adata_extract_chunk, dataset_name=dataset_name)


def test_extract_filtered_by_homo_sapiens():
    """
    Tests data extraction for Homo sapiens. Initiates extraction tables, extracts data, verifies chunks, and
    cleans resources. Focuses on the Homo sapiens organism.

    :raises AssertionError: If extracted chunks are not verified.
    """
    test_extract_data_dir = f"{constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX}__data"
    num_extract_chunks_to_check = 9
    logger.info("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Homo sapiens",
    )
    logger.info("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=num_extract_chunks_to_check - 1,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
        obs_columns_to_include_str=constants.OBS_COLUMNS_TO_INCLUDE,
    )
    logger.info("Verifying extract files...")
    gcs_extract_file_paths = [
        f"{test_extract_data_dir}/{constants.GCS_EXTRACT_FILE_NAME_PREFIX.format(chunk_id=x)}"
        for x in range(num_extract_chunks_to_check)
    ]
    verify_extracted_chunks(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_file_paths=gcs_extract_file_paths,
        dataset_name=constants.DATASET_NAME,
    )
    logger.info("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_mus_mus():
    """
    Tests data extraction for Mus musculus. Initiates extraction tables, extracts data, verifies chunks, and
    cleans resources. Focuses on the Mus musculus organism.

    :raises AssertionError: If extracted chunks are not verified.
    """
    test_extract_data_dir = f"{constants.MUS_MUS_EXTRACT_TABLE_PREFIX}__data"
    num_extract_chunks_to_check = 3
    logger.info("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.MUS_MUS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Mus musculus",
    )
    logger.info("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=num_extract_chunks_to_check - 1,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
        obs_columns_to_include_str=constants.OBS_COLUMNS_TO_INCLUDE,
    )
    logger.info("Verifying extract files...")
    gcs_extract_file_paths = [
        f"{test_extract_data_dir}/{constants.GCS_EXTRACT_FILE_NAME_PREFIX.format(chunk_id=x)}"
        for x in range(num_extract_chunks_to_check)
    ]
    verify_extracted_chunks(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_MUS_MUS,
        gcs_extract_file_paths=gcs_extract_file_paths,
        dataset_name=constants.DATASET_NAME,
    )
    logger.info("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        extract_table_prefix=constants.MUS_MUS_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_homo_sapiens_small_chunk_size():
    """
    Test data extraction for Homo sapiens with smaller chunk size. This function prepares extraction tables, extracts
    data, verifies the resulting chunks, and then cleans up the resources. It focuses on data for Homo sapiens using a
    smaller chunk size for extraction.

    :raises AssertionError: If extracted chunks are not verified.
    """
    test_extract_data_dir = f"{constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX}__data"
    num_extract_chunks_to_check = 18
    logger.info("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=5000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_organism="Homo sapiens",
    )
    logger.info("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=num_extract_chunks_to_check - 1,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
        obs_columns_to_include_str=constants.OBS_COLUMNS_TO_INCLUDE,
    )
    logger.info("Verifying extract files...")
    gcs_extract_file_paths = [
        f"{test_extract_data_dir}/{constants.GCS_EXTRACT_FILE_NAME_PREFIX.format(chunk_id=x)}"
        for x in range(num_extract_chunks_to_check)
    ]
    verify_extracted_chunks(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_file_paths=gcs_extract_file_paths,
        dataset_name=constants.DATASET_NAME,
    )
    logger.info("Cleaning up infrastructure from files that were produced by the test...")
    clean_up_cloud_from_test_case(
        bucket_name=constants.GCS_BUCKET_NAME,
        extract_table_prefix=constants.HOMO_SAPIENS_5k_EXTRACT_TABLE_PREFIX,
    )


def test_extract_filtered_by_datasets():
    """
    Test data extraction filtered by specific datasets. This function prepares extraction tables, extracts data based
    on specific datasets, verifies the resulting chunks, and cleans up afterwards. It emphasizes filtering based on
    datasets.

    :raises AssertionError: If extracted chunks are not verified.
    """
    test_extract_data_dir = f"{constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX}__data"
    num_extract_chunks_to_check = 6
    logger.info("Preparing extract tables...")
    prepare_extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX,
        min_observed_cells=0,
        fq_allowed_original_feature_ids=constants.HOMO_SAPIENS_GENE_SCHEMA,
        extract_bin_size=10000,
        bucket_name=constants.GCS_BUCKET_NAME,
        filter_by_datasets="tg66rtgh7y-test-source-file-1_ingest_info.avro,yt55tgy4o-test-source-file-3_ingest_info.avro",
    )
    logger.info("Extracting data...")
    extract(
        dataset=constants.DATASET_NAME,
        extract_table_prefix=constants.FILTER_BY_DATASET_EXTRACT_TABLE_PREFIX,
        start_bin=0,
        end_bin=num_extract_chunks_to_check - 1,
        output_bucket_name=constants.GCS_BUCKET_NAME,
        output_bucket_directory=test_extract_data_dir,
        obs_columns_to_include_str=constants.OBS_COLUMNS_TO_INCLUDE,
    )
    logger.info("Verifying extract files...")
    gcs_extract_file_paths = [
        f"{test_extract_data_dir}/{constants.GCS_EXTRACT_FILE_NAME_PREFIX.format(chunk_id=x)}"
        for x in range(num_extract_chunks_to_check)
    ]
    verify_extracted_chunks(
        gcs_bucket_name=constants.GCS_BUCKET_NAME,
        gcs_input_file_paths=constants.GCS_INPUT_PATHS_HOMO_SAPIENS,
        gcs_extract_file_paths=gcs_extract_file_paths,
        dataset_name=constants.DATASET_NAME,
    )
