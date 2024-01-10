<%
    from casp.workflows.kubeflow.machine_specs import TDigestMachineSpec

    accelerator = "cpu" if TDigestMachineSpec.accelerator_type is None else "gpu"
    num_nodes = TDigestMachineSpec.replica_count
    devices = 4 if TDigestMachineSpec.accelerator_count is None else TDigestMachineSpec.accelerator_count
%>
# lightning.pytorch==2.0.9
seed_everything: true
trainer:
  accelerator: cpu
  callbacks:
    - class_path: lightning.pytorch.callbacks.ModelCheckpoint
      init_args:
        filename: 'model'
  strategy:
    class_path: lightning.pytorch.strategies.DDPStrategy
    init_args:
      timeout: 0:30:00
      start_method: popen
      dim: 0
      broadcast_buffers: false
      bucket_cap_mb: 25
      find_unused_parameters: false
      check_reduction: false
      gradient_as_bucket_view: false
      static_graph: false
  devices: ${devices}
  num_nodes: ${num_nodes}
  precision: 32-true
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
  default_root_dir: /gcs/cellarium-file-system/curriculum/${curriculum_name}/models/${model_name}

model:
  model: cellarium.ml.models.TDigest
  transforms:
    - class_path: cellarium.ml.transforms.NormalizeTotal
      init_args:
        target_count: ${normalize_total_target_count}
        eps: ${normalize_total_eps}
% if use_log_1p:
    - class_path: cellarium.ml.transforms.Log1p
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