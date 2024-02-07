<%
    from casp.workflows.kubeflow.machine_specs import PCAEmbedMachineSpec

    accelerator = "cpu" if PCAEmbedMachineSpec.accelerator_type is None else "gpu"
    num_nodes = PCAEmbedMachineSpec.replica_count
    devices = 4 if PCAEmbedMachineSpec.accelerator_count is None else PCAEmbedMachineSpec.accelerator_count
%>
# lightning.pytorch==2.0.9
seed_everything: true
trainer:
  accelerator: gpu
  strategy: auto
  devices: ${devices}
  num_nodes: ${num_nodes}
  precision: 32-true
  callbacks:
  - class_path: cellarium.ml.callbacks.PredictionWriter
    init_args:
      output_dir: /gcs/cellarium-file-system/curriculum/${curriculum_name}/models/${model_name}/embeddings
      prediction_size: ${embedding_dimension}
  fast_dev_run: false
  max_epochs: 1
  max_steps: -1
  overfit_batches: 0.0
  check_val_every_n_epoch: 1
  accumulate_grad_batches: 1
  inference_mode: true
  use_distributed_sampler: true
  detect_anomaly: false
  barebones: false
  sync_batchnorm: false
  reload_dataloaders_every_n_epochs: 0

model:
  model:
    class_path: cellarium.ml.models.IncrementalPCA
    init_args:
      n_components: ${n_components}
      svd_lowrank_niter: 2
      perform_mean_correction: false
  transforms:
      - class_path: cellarium.ml.transforms.NormalizeTotal
        init_args:
          target_count: ${normalize_total_target_count}
          eps: ${normalize_total_eps}
% if use_log_1p:
      - class_path: cellarium.ml.transforms.Log1p
% endif
% if use_filter:
      - class_path: cellarium.ml.transforms.Filter
        init_args:
          filter_list:
            !FileLoader
            file_path: ${filter_file_path}
            loader_fn: pandas.read_csv
            attr: original_feature_id
            convert_fn: pandas.Series.to_list
% endif
data:
  dadc:
    class_path: cellarium.ml.data.DistributedAnnDataCollection
    init_args:
      filenames: gs://cellarium-file-system/curriculum/${curriculum_name}/extract_files/extract_{${extract_start}..${extract_end}}.h5ad
      shard_size: ${shard_size}
      last_shard_size: ${last_shard_size}
      max_cache_size: 2
      cache_size_strictly_enforced: true
      indices_strict: true
      obs_columns_to_validate: ["total_mrna_umis"]
  batch_size: 10_000
  shuffle: false
  seed: 0
  drop_last: true
  test_mode: false
  num_workers: 4
  batch_keys:
    x_ng:
      attr: X
      convert_fn: cellarium.ml.utilities.data.densify
    var_names_g:
      attr: var_names
    total_mrna_umis_n:
      attr: obs
      key: total_mrna_umis
return_predictions: false
ckpt_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${model_name}/lightning_logs/version_0/checkpoints/model.ckpt