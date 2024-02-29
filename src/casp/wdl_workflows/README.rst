WDL Workflows
=============

This module contains all the instructions for automated execution of cas services.

Module Structure
----------------
- bq_ops
   - anndata_to_ingest_files
   - ingest_files_to_bq
   - prepare_extract
   - extract
   - test_data_cycle
- data_embedding

Set Up
------

Cromwell Execution Manager
~~~~~~~~~~~~~~~~~~~~~~~~~~
To execute these scripts in Cloud Environment you would need a cromwell task execution manager.
Broad already has a cromwell server running in the cloud. You can use that server to execute these scripts.
To interact with this server you would need to install `cromshell` (cromwell cli tool) on your local machine.
You can do that by this command:

.. code:: bash

    $ pip install cromshell

After installation, you would need to specify a cromwell server URL. You can do that by this command:

.. code:: bash

    $ cromshell config set url ${CROMWELL_SERVER_URL}

Please refer to the `Cromshell Github <https://github.com/broadinstitute/cromshell>`__ page for more

VPN
~~~
To interact with the cromwell server you would need to be on the Broad
VPN. Make you sure that you’re using ``Z-Duo-Broad-NonSplit-VPN`` VPN mode, otherwise cromwell server won’t be accessible. Please refer to the
`Broad VPN <https://intranet.broadinstitute.org/bits/service-catalog/networking/vpn>`__ page for more information.

Usage
~~~~~
To execute a workflow you would need to run this command:

.. code:: bash

    $ cromshell submit ${WDL_FILE_PATH} ${INPUT_JSON_FILE_PATH}


To check the status of your workflow runs you can run this command:

.. code:: bash

    $ cromshell list -u -c

Workflow Inputs Files
---------------------

More detailed description of what these inputs are, please look at the `casp.services.bq_ops` directory and its subdirectories for each individual operation.

bq_ops.anndata_to_ingest_files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters
...........
- **CASAnndataToIngestFiles.anndata_to_ingest_files.docker_image** - Docker image to use for this operation.
- **CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_stage_dir** - GCS directory to stage the output files.
- **CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_input_bucket** - Working GCS Bucket name.
- **CASAnndataToIngestFiles.anndata_to_ingest_files.original_feature_id_lookup** - A column name in var dataframe from where to get original feature ids. In most of the cases it will be a column with ENSEMBL gene IDs. if `index`, then the index of the dataframe will be used.
- **CASAnndataToIngestFiles.convert_args** - List of dictionaries with the following keys |br|
  **df_filename** - A path to the anndata file in GCS bucket |br|
  **cas_cell_index** - A starting index for the cells in the output files |br|
  **cas_feature_index** - A starting index for the features in the output files |br|

Example
.......

.. code:: json

    {
       "CASAnndataToIngestFiles.anndata_to_ingest_files.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
       "CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_stage_dir": "cromwell_50m",
       "CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_input_bucket": "cellarium-file-system",
       "CASAnndataToIngestFiles.anndata_to_ingest_files.original_feature_id_lookup": "index",
       "CASAnndataToIngestFiles.convert_args": [
          {
             "df_filename": "census_data/0129dbd9-a7d3-4f6b-96b9-1da155a93748-census-dataset.h5ad",
             "cas_cell_index": 0,
             "cas_feature_index": 0
          },
          {
             "df_filename": "census_data/04a23820-ffa8-4be5-9f65-64db15631d1e-census-dataset.h5ad",
             "cas_cell_index": 1000000,
             "cas_feature_index": 1000000
          }
       ]
    }

bq_ops.ingest_file_to_bq
~~~~~~~~~~~~~~~~~~~~~~~~

Parameters
...........
- **CASIngestFilesToBQ.ingest_files_to_bq.docker_image** - Docker image to use for this operation.
- **CASIngestFilesToBQ.ingest_files_to_bq.gcs_stage_dir** - GCS directory to stage the output files.
- **CASIngestFilesToBQ.ingest_files_to_bq.gcs_bucket_name** - Working GCS Bucket name.
- **CASIngestFilesToBQ.ingest_files_to_bq.dataset** - BigQuery dataset name where to ingest the data. If dataset doesn't exist, it will be created.

Example
.......

.. code:: json

    {
      "CASIngestFilesToBQ.ingest_files_to_bq.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
      "CASIngestFilesToBQ.ingest_files_to_bq.gcs_stage_dir": "cromwell_test_10k",
      "CASIngestFilesToBQ.ingest_files_to_bq.gcs_bucket_name": "cellarium-file-system",
      "CASIngestFilesToBQ.ingest_files_to_bq.dataset": "cas_test_dataset"
    }


bq_ops.precalculate_fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters
..........
- **CASPrecalculateFields.precalculate_fields.docker_image** - Docker image to use for this operation
- **CASPrecalculateFields.precalculate_fields.dataset** - BigQuery dataset name where to ingest the data.
- **CASPrecalculateFields.precalculate_fields.fields** - A comma separated list of fields to precalculate. Currently only `total_mrna_umis` is supported.

Example
.......

.. code:: json

    {
       "CASPrecalculateFields.precalculate_fields.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
       "CASPrecalculateFields.precalculate_fields.dataset": "cas_test_dataset",
       "CASPrecalculateFields.precalculate_fields.fields": "total_mrna_umis"
    }


bq_ops.prepare_extract
~~~~~~~~~~~~~~~~~~~~~~
Current input file could have a filed `CASPrepareExtractBQ.prepare_extract.filters_json_path`. This is a gs json path. Please find an example of this filter file attached as well

Parameters
...........
- **CASPrepareExtractBQ.prepare_extract.docker_image** - Docker image to use for this operation.
- **CASPrepareExtractBQ.prepare_extract.bq_dataset** - BigQuery dataset name where to ingest the data.
- **CASPrepareExtractBQ.prepare_extract.extract_table_prefix** - Prefix for the extract table name. The final extract tables name will be named like `${extract_table_prefix}_cas_cell_info`
- **CASPrepareExtractBQ.prepare_extract.extract_bin_size** - Size of the bin for the extract table, usually we put 10000
- **CASPrepareExtractBQ.prepare_extract.bucket_name** - Working GCS Bucket name
- **CASPrepareExtractBQ.prepare_extract.obs_columns_to_include** - A comma separated list of columns to include in the extract table.
- **CASPrepareExtractBQ.prepare_extract.fq_allowed_original_feature_ids** - A fully qualified table name of the reference data table with the feature schema.
- **CASPrepareExtractBQ.prepare_extract.extract_bucket_path** - A GCS path to the extract table. This path is used for creating metadata files for the extract script.
- **CASPrepareExtractBQ.prepare_extract.filters_json_path** - A GCS path to a json file with filters. Please find an example of this filter file attached as well

Example
.......

.. code:: json

    {
      "CASPrepareExtractBQ.prepare_extract.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
      "CASPrepareExtractBQ.prepare_extract.bq_dataset": "cas_test_dataset",
      "CASPrepareExtractBQ.prepare_extract.extract_table_prefix": "fg_extract",
      "CASPrepareExtractBQ.prepare_extract.extract_bin_size": 10000,
      "CASPrepareExtractBQ.prepare_extract.bucket_name": "cellarium-file-system",
      "CASPrepareExtractBQ.prepare_extract.extract_bucket_path": "curriculum/fg_extract",
      "CASPrepareExtractBQ.prepare_extract.filters_json_path": "gs://cellarium-file-system/curriculum/extract_filters/filters_mus_mus_brain.json",
      "CASPrepareExtractBQ.prepare_extract.fq_allowed_original_feature_ids": "dsp-cell-annotation-service.cas_reference_data.refdata-gex-GRCh38-2020-A",
      "CASPrepareExtractBQ.prepare_extract.obs_columns_to_include": "cell_type,total_mrna_umis,donor_id,assay,development_stage,disease,organism,sex,tissue"
    }

An example of JSON object for `CASPrepareExtractBQ.prepare_extract.filters_json_path` (you'd need to put this file in a GCS bucket and provide a path to the workflow input file):

.. code:: json

    {
      "organism__eq": "Mus musculus",
      "cell_type__in": ["L6b glutamatergic cortical neuron", "interneuron", "inhibitory interneuron", "cerebellar Golgi cell"],
      "is_primary_data__eq": true
    }

.. note::
    Constructing filters

    Filters is a dictionary containing filter criteria, structured as ``{"column_name__filter_type": "value"}``.

    Supported filter_types |br|
    ``"eq"`` - Used for an 'equals' comparison. |br|
    `Example`: ``{"organism__eq": "Homo sapiens"}`` results in ``organism='Homo sapiens'``.

    ``"in"`` - Used for an 'in' comparison with a set of values. |br|
    `Example`: ``{"cell_type__in": ["T cell", "neuron"]}`` results in ``cell_type in ('T cell', 'neuron')``.

bq_ops.extract
~~~~~~~~~~~~~~

Parameters
..........

- **CASExtractBQ.extract.docker_image** - Docker image to use for this operation.
- **CASExtractBQ.extract.bq_dataset** - BigQuery dataset name where to ingest the data.
- **CASExtractBQ.extract.extract_table_prefix** - Prefix for the extract table name. This value should be the same as the one used in `CASPrepareExtractBQ.prepare_extract.extract_table_prefix`
- **CASExtractBQ.bin_borders** - A list of lists with bin borders. Each list should contain two numbers, the first one is the start of the bin, the second one is the end of the bin.

For example, `[[0, 9], [10, 19], [20, 29]]` will create 30 bins: 0-9, 10-19, 20-29. Each bin will be a
separate `.h5ad` extract file. Each bin group will be executed in parallel on a separate VM. Number of bins per group
should correspond to a number of CPU cores in each of the machine. It is not recommended to have more than 50 groups at
the same time because cromwell wouldn't be able to manage all the machines and will lose some of the groups.
- **CASExtractBQ.extract.output_bucket_name** - Working GCS Bucket name
- **CASExtractBQ.extract.extract_bucket_path** - A GCS path to the extract table. Should be the same as the one used in `CASPrepareExtractBQ.prepare_extract.extract_bucket_path` to use metadata file produced by `extract` script


Example
.......

.. code:: json

    {
      "CASExtractBQ.extract.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
      "CASExtractBQ.extract.bq_dataset": "cas_test_dataset",
      "CASExtractBQ.extract.extract_table_prefix": "fg_extract",
      "CASExtractBQ.bin_borders": [[0, 9], [10, 19], [20, 29], [30, 39], [40, 49], [50, 59]],
      "CASExtractBQ.extract.output_bucket_name": "cellarium-file-system",
      "CASExtractBQ.extract.extract_bucket_path": "curriculum/fg_extract",
      "CASExtractBQ.extract.obs_columns_to_include": "cell_type,total_mrna_umis,donor_id,assay,development_stage,disease,organism,sex,tissue"
    }


bq_ops.test_data_cycle
~~~~~~~~~~~~~~~~~~~~~~
Only requires a docker image

.. code:: json

    {
       "CASTestDataCycle.test_data_cycle.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1"
    }

model_training.train_incremental_pca (Deprecated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: json

    {
      "CASTrainIncrementalPCA.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
      "CASTrainIncrementalPCA.bucket_name": "dsp-cell-annotation-service",
      "CASTrainIncrementalPCA.data_storage_path": "cas_50m_homo_sapiens_extract_4m",
      "CASTrainIncrementalPCA.checkpoint_save_path": "pca_incremental_4m_june",
      "CASTrainIncrementalPCA.n_components": 512,
      "CASTrainIncrementalPCA.batch_size": 10000,
      "CASTrainIncrementalPCA.use_gpu": true
    }

data_embedding (Deprecated)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: json

    {
      "CASPEmbedData.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
      "CASPEmbedData.bucket_name": "dsp-cell-annotation-service",
      "CASPEmbedData.data_storage_path": "cas_50m_homo_sapiens_extract_4m",
      "CASPEmbedData.dm_storage_path": "models/incremental_pca_003.pickle",
      "CASPEmbedData.output_storage_path": "embeddings_incremental_pca_003",
      "CASPEmbedData.running_script": "casp/services/data_embedding/main.py"
    }