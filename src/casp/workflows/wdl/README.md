# CAS WDL Workflows Directory
This module contains all the instructions for automated execution of cas services.
## Module Structure
 - bq_ops
   - anndata_to_ingest_files
   - ingest_files_to_bq
   - prepare_extract
   - extract
   - test_data_cycle
 - data_embedding
 - model_training
   - train_incremental_pca

### Workflow Inputs
More detailed description of what these inputs are, please look at the `casp.services` directory and its subdirectories for each individual operation
### Input Examples
bq_ops.anndata_to_ingest_files:
```JSON
{
   "CASAnndataToIngestFiles.anndata_to_ingest_files.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
   "CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_stage_dir": "cromwell_50m",
   "CASAnndataToIngestFiles.anndata_to_ingest_files.gcs_input_bucket": "dsp-cell-annotation-service",
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
 ```
bq_ops.ingest_file_to_bq:
```JSON
{
  "CASIngestFilesToBQ.ingest_files_to_bq.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
  "CASIngestFilesToBQ.ingest_files_to_bq.gcs_stage_dir": "cromwell_test_10k",
  "CASIngestFilesToBQ.ingest_files_to_bq.gcs_bucket_name": "dsp-cell-annotation-service",
  "CASIngestFilesToBQ.ingest_files_to_bq.dataset": "cas_test_dataset"
}
```
bq_ops.prepare_extract:
Current input file could have a filed `CASPrepareExtractBQ.prepare_extract.filters_json_path`. This is a gs json path. Please find an example of this filter file attached as well 
```JSON
{
  "CASPrepareExtractBQ.prepare_extract.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
  "CASPrepareExtractBQ.prepare_extract.bq_dataset": "cas_test_dataset",
  "CASPrepareExtractBQ.prepare_extract.extract_table_prefix": "fg_extract",
  "CASPrepareExtractBQ.prepare_extract.extract_bin_size": 10000,
  "CASPrepareExtractBQ.prepare_extract.bucket_name": "dsp-cell-annotation-service",
  "CASPrepareExtractBQ.prepare_extract.filters_json_path": "gs://cellarium-file-system/curriculum/extract_filters/filters_mus_mus_brain.json",
  "CASPrepareExtractBQ.prepare_extract.obs_columns_to_include_str": "cell_type,total_mrna_umis,donor_id,assay,development_stage,disease,organism,sex,tissue"
}
```
CASPrepareExtractBQ.prepare_extract.filters_json_path example:
```JSON
{
  "organism__eq": "Mus musculus",
  "cell_type__in": ["L6b glutamatergic cortical neuron", "interneuron", "inhibitory interneuron", "cerebellar Golgi cell"],
  "is_primary_data__eq": true
}
```
bq_ops.extract:
```JSON
{
  "CASExtractBQ.extract.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
  "CASExtractBQ.extract.project_id": "dsp-cell-annotation-service",
  "CASExtractBQ.extract.bq_dataset": "cas_test_dataset",
  "CASExtractBQ.extract.extract_table_prefix": "fg_extract",
  "CASExtractBQ.bin_borders": [[0, 9], [10, 19], [20, 29], [30, 39], [40, 49], [50, 59]],
  "CASExtractBQ.extract.output_bucket_name": "dsp-cell-annotation-service",
  "CASExtractBQ.extract.output_bucket_directory": "cas_test_dataset_extract", 
  "CASExtractBQ.extract.obs_columns_to_include_str": "cell_type,total_mrna_umis,donor_id,assay,development_stage,disease,organism,sex,tissue"
}
```
bq_ops.test_data_cycle (requires only docker image):
```JSON
{
   "CASTestDataCycle.test_data_cycle.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1"
}
```
bq_ops.precalculate_fields 
```JSON
{
   "CASPrecalculateFields.precalculate_fields.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
   "CASPrecalculateFields.precalculate_fields.dataset": "cas_test4m_dataset",
   "CASPrecalculateFields.precalculate_fields.fields": "total_mrna_umis"
}
```
model_training.train_incremental_pca:
```JSON
{
  "CASTrainIncrementalPCA.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
  "CASTrainIncrementalPCA.bucket_name": "dsp-cell-annotation-service",
  "CASTrainIncrementalPCA.data_storage_path": "cas_50m_homo_sapiens_extract_4m",
  "CASTrainIncrementalPCA.checkpoint_save_path": "pca_incremental_4m_june",
  "CASTrainIncrementalPCA.n_components": 512,
  "CASTrainIncrementalPCA.batch_size": 10000,
  "CASTrainIncrementalPCA.use_gpu": true
}
```
data_embedding:
```JSON
{
  "CASPEmbedData.docker_image": "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:1.0a1",
  "CASPEmbedData.bucket_name": "dsp-cell-annotation-service",
  "CASPEmbedData.data_storage_path": "cas_50m_homo_sapiens_extract_4m",
  "CASPEmbedData.dm_storage_path": "models/incremental_pca_003.pickle",
  "CASPEmbedData.output_storage_path": "embeddings_incremental_pca_003",
  "CASPEmbedData.running_script": "casp/services/data_embedding/main.py"
}
```