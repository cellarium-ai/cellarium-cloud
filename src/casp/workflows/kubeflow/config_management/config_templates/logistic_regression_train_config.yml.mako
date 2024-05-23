<%
    from casp.workflows.kubeflow import machine_specs_utils

    machine_specs_info = machine_specs_utils.read_machine_specs()
    current_spec_info = machine_specs_info["logistic_regression_train"]
    accelerator = "cpu" if current_spec_info["accelerator_type"] is None else "gpu"
    num_nodes = current_spec_info["replica_count"]
    devices = 4 if current_spec_info["accelerator_count"] is None else current_spec_info["accelerator_count"]
%>
# lightning.pytorch==2.0.9
seed_everything: true
trainer:
  accelerator: ${accelerator}
  strategy: auto
  devices: ${devices}
  num_nodes: ${num_nodes}
  precision: 32-true
  fast_dev_run: false
  max_epochs: ${max_epochs}
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
  default_root_dir: /gcs/cellarium-file-system/curriculum/${curriculum_name}/models/${model_name}
  callbacks:
  - class_path: lightning.pytorch.callbacks.LearningRateMonitor
    init_args:
      logging_interval: step
      log_momentum: false
  - class_path: lightning.pytorch.callbacks.ModelCheckpoint
    init_args:
      every_n_train_steps: 1000
      save_top_k: -1

model:
  model:
    class_path: cellarium.ml.models.LogisticRegression
    init_args:
      W_prior_scale: ${W_prior_scale}
      W_init_scale: 1.0
      seed: 0
  optim_fn: torch.optim.Adam
  optim_kwargs:
    lr: ${lr}
  scheduler_fn: cellarium.ml.lr_schedulers.LinearLR
  scheduler_kwargs:
    num_warmup_steps: 0
    num_training_steps: ${scheduler_num_training_steps}
  transforms:
    - class_path: cellarium.ml.transforms.NormalizeTotal
      init_args:
        target_count: ${normalize_total_target_count}
        eps: ${normalize_total_eps}
% if use_log_1p:
    - class_path: cellarium.ml.transforms.Log1p
% endif
% if use_z_score:
    - class_path: cellarium.ml.transforms.ZScore
      init_args:
        mean_g:
          !CheckpointLoader
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${mean_var_std_model_name}/lightning_logs/version_0/checkpoints/model.ckpt
          attr: model.mean_g
          convert_fn: null
        std_g:
          !CheckpointLoader
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${mean_var_std_model_name}/lightning_logs/version_0/checkpoints/model.ckpt
          attr: model.std_g
          convert_fn: null
        var_names_g:
          !CheckpointLoader
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${mean_var_std_model_name}/lightning_logs/version_0/checkpoints/model.ckpt
          attr: model.var_names_g
          convert_fn: numpy.ndarray.tolist
% endif
% if use_divide_by_scale:
    - class_path: cellarium.ml.transforms.DivideByScale
      init_args:
        scale_g:
          !CheckpointLoader
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${tdigest_model_name}/lightning_logs/version_0/checkpoints/model.ckpt
          attr: model.median_g
          convert_fn: torch.Tensor.float
        var_names_g:
          !CheckpointLoader
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/${tdigest_model_name}/lightning_logs/version_0/checkpoints/model.ckpt
          attr: model.var_names_g
          convert_fn: numpy.ndarray.tolist
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
      obs_columns_to_validate: ["total_mrna_umis", "cell_type"]
  batch_size: ${batch_size}
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
    y_n:
      attr: obs
      key: cell_type
      convert_fn: cellarium.ml.utilities.data.categories_to_codes
    total_mrna_umis_n:
      attr: obs
      key: total_mrna_umis