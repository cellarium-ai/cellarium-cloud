# BQ Ops Services
=================

These are wrappers around the code in `casp.bq_scripts`. These wrappers try to follow the conventions of the old scripts from `casp.bq_scripts` to avoid massive code base change in one time. Each of the `bq_ops` services has its own Cromwell workflow in `casp.cromwell.bq_ops`.

Prerequisites
-------------

1. Google BigQuery API Allowed
2. Credentials dumped into `src/casp/services/.env`
3. `cas_reference_data` dataset in BigQuery with tables for each schema

Running Locally
---------------

Prepare Extract
~~~~~~~~~~~~~~~


.. list-table:: Environment Variables
   :header-rows: 1

   * - Variable Name
     - Example
     - Description
   * - dataset_name
     - your_dataset_name
     - BigQuery dataset
   * - extract_table_prefix
     - extract_table_prefix_name
     - Name used by the script to aggregate the data
   * - bucket_name
     - bucket-name
     - Working bucket name
   * - extract_bin_size
     - 10000
     - Size of extract bin (anndata file batch size)
   * - fq_allowed_original_feature_ids
     - fqn
     - Fully qualified name of the table with a list of Ensembl gene IDs
   * - filters_json_path
     - ./filters.json
     - Path to a JSON file with filters (can be local or gs:// path)
   * - extract_bucket_path
     - path/to/dir
     - Path on the bucket where this extract should be stored
   * - obs_columns_to_include
     - cell_type,cell_type_ontology_term_id
     - Comma-separated list of column names to include in obs from BigQuery cell info



.. code-block:: bash
    $ python casp/services/bq_ops/prepare_extract/main.py \
    --dataset="cas_50m_dataset" \
    --extract_table_prefix="mus_mus_test_may_14" \
    --bucket_name="cellarium-file-system" \
    --extract_bin_size=10000 \
    --fq_allowed_original_feature_ids="dsp-cell-annotation-service.cas_reference_data.refdata-gex-mm10-2020-A" \
    --filters_json_path="./filters.json" \
    --extract_bucket_path="curriculum/curriculum_name" \
    --obs_columns_to_include="cell_type,cell_type_ontology_term_id,disease,disease_ontology_term_id,total_mrna_umis,i.dataset_id"



Test Data Cycle
~~~~~~~~~~~~~~~

A script that tests data infrastructure if it could produce relevant chunks after harmonization and randomization. It is excluded from normal tests as it is not a unit test. Just run this test manually:

.. code-block:: bash

    $ python casp/services/bq_ops/test_data_cycle/main.py

You can run this as a WDL workflow. The test raises an error after all test cases are finished. This would ensure your WDL workflow would become marked as `Failed`.
