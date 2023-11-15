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
      prediction_size: ${prediction_size}
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
    class_path: cellarium.ml.models.IncrementalPCAFromCLI
    init_args:
      k_components: ${k_components}
      svd_lowrank_niter: 2
      perform_mean_correction: false
      target_count: 10000
  default_lr: 0.001
data:
  filenames: gs://cellarium-file-system/curriculum/${curriculum_name}/extract_files/extract_{${extract_start}..${extract_end}}.h5ad
  shard_size: ${shard_size}
  last_shard_size: ${last_shard_size}
  max_cache_size: 2
  cache_size_strictly_enforced: true
  indices_strict: true
  batch_size: 10000
  shuffle: false
  seed: 0
  drop_last: true
  test_mode: false
  num_workers: 4
  batch_keys:
    X:
      attr: X
      convert_fn: cellarium.ml.utilities.data.densify
    var_names:
      attr: var_names
    obs_names:
      attr: obs_names
return_predictions: false
ckpt_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${model_name}/lightning_logs/version_0/checkpoints/model.ckpt