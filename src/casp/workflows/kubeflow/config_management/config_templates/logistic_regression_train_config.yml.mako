<%
    from casp.workflows.kubeflow.machine_specs import LogisticRegressionTrainMachineSpec

    accelerator = "cpu" if LogisticRegressionTrainMachineSpec.accelerator_type is None else "gpu"
    num_nodes = LogisticRegressionTrainMachineSpec.replica_count
    devices = 4 if LogisticRegressionTrainMachineSpec.accelerator_count is None else LogisticRegressionTrainMachineSpec.accelerator_count
%>
# lightning.pytorch==2.0.9
seed_everything: true
trainer:
  accelerator: ${accelerator}
  callbacks:
    - class_path: lightning.pytorch.callbacks.ModelCheckpoint
      init_args:
        filename: 'model'
  strategy: auto
  logger:
    class_path: lightning.pytorch.loggers.NeptuneLogger
    init_args:
      api_key: eyJhcGlfYWRkcmVzcyI6Imh0dHBzOi8vYXBwLm5lcHR1bmUuYWkiLCJhcGlfdXJsIjoiaHR0cHM6Ly9hcHAubmVwdHVuZS5haSIsImFwaV9rZXkiOiJkNjNhOTY3ZS1mOGU2LTQ2ZGItYTFmOS01MGY4ZDdiNGU1YTcifQ==
      project: 'cellarium-pipelines/cas'
      name: '${curriculum_name}-${model_name}'

  devices: ${devices}
  num_nodes: ${num_nodes}
  precision: 32-true
  fast_dev_run: false
  max_epochs: 3
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

model:
  model:
    class_path: cellarium.ml.models.LogisticRegression
    init_args:
          W_prior_scale: 1.0
          W_init_scale: 1.0
          seed: 0
  optim_fn: torch.optim.Adam
  optim_kwargs:
    lr: 0.001
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
          file_path: gs://cellarium-file-system/curriculum/${curriculum_name}/models/shared_meta/${tdigest_model_name}-non-nan-features.csv
          loader_fn: pandas.read_csv
          attr: original_feature_id
          convert_fn: pandas.Series.to_list
% endif

data:
  filenames: gs://cellarium-file-system/curriculum/${curriculum_name}/extract_files/extract_{${extract_start}..${extract_end}}.h5ad
  shard_size: ${shard_size}
  last_shard_size: ${last_shard_size}
  max_cache_size: 2
  cache_size_strictly_enforced: true
  indices_strict: true
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
    y_n:
      attr: obs
      key: cell_type
      convert_fn: cellarium.ml.utilities.data.categories_to_codes
    total_mrna_umis_n:
      attr: obs
      key: total_mrna_umis