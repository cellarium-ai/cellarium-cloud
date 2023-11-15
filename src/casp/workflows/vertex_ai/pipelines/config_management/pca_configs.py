import typing as t

import yaml
from mako.template import Template
from smart_open import open

from casp.services import settings

PIPELINE_CONFIGS_TEMPLATE_PATH = f"{settings.CAS_DIR}/workflows/vertex_ai/pipelines/config_management/config_templates"

TRAIN_TEMPLATE_PATH = f"{PIPELINE_CONFIGS_TEMPLATE_PATH}/ipca_train_config.yml.mako"
EMBED_TEMPLATE_PATH = f"{PIPELINE_CONFIGS_TEMPLATE_PATH}/ipca_embed_config.yml.mako"
REGISTER_MODEL_TEMPLATE_PATH = f"{PIPELINE_CONFIGS_TEMPLATE_PATH}/ipca_register_model_config.yml.mako"

CONFIG_BUCKET_PATH = "gs://cellarium-file-system/ml-configs/ipca"

TRAIN_EMBED_REQUIRED_KEYS = [
    "curriculum_name",
    "model_name",
    "k_components",
    "shard_size",
    "train_devices",
    "train_num_nodes",
    "train_extract_start",
    "train_extract_end",
    "train_last_shard_size",
    "embed_devices",
    "embed_num_nodes",
    "embed_extract_start",
    "embed_extract_end",
    "embed_last_shard_size",
    "embed_prediction_size",
]


def _validate_train_embed_config(train_embed_config: t.Dict[str, t.Any]):
    for key in TRAIN_EMBED_REQUIRED_KEYS:
        if key not in train_embed_config:
            raise KeyError(f"Missing required key {key} in train_embed_config")


def create_train_embed_configs(
    curriculum_name: str,
    model_name: str,
    shard_size: int,
    k_components: int,
    train_devices: int,
    train_num_nodes: int,
    train_extract_start: int,
    train_extract_end: int,
    train_last_shard_size: int,
    embed_devices: int,
    embed_num_nodes: int,
    embed_extract_start: int,
    embed_extract_end: int,
    embed_last_shard_size: int,
    embed_prediction_size: int,
) -> t.Tuple[str, str, str]:
    train_config_template = Template(filename=TRAIN_TEMPLATE_PATH)
    embed_config_template = Template(filename=EMBED_TEMPLATE_PATH)
    register_model_config_template = Template(filename=REGISTER_MODEL_TEMPLATE_PATH)

    rendered_train_config = train_config_template.render(
        curriculum_name=curriculum_name,
        model_name=model_name,
        shard_size=shard_size,
        devices=train_devices,
        num_nodes=train_num_nodes,
        extract_start=train_extract_start,
        extract_end=train_extract_end,
        last_shard_size=train_last_shard_size,
        k_components=k_components,
    )
    rendered_embed_config = embed_config_template.render(
        curriculum_name=curriculum_name,
        model_name=model_name,
        shard_size=shard_size,
        devices=embed_devices,
        num_nodes=embed_num_nodes,
        extract_start=embed_extract_start,
        extract_end=embed_extract_end,
        last_shard_size=embed_last_shard_size,
        k_components=k_components,
        prediction_size=embed_prediction_size,
    )
    rendered_register_config = register_model_config_template.render(
        curriculum_name=curriculum_name, model_name=model_name, embedding_dimension=embed_prediction_size
    )
    train_config_path = f"{CONFIG_BUCKET_PATH}/{curriculum_name}-{model_name}-train-config.yml"
    embed_config_path = f"{CONFIG_BUCKET_PATH}/{curriculum_name}-{model_name}-embed-config.yml"
    register_config_path = f"{CONFIG_BUCKET_PATH}/{curriculum_name}-{model_name}-register-config.yml"

    with open(train_config_path, "w") as f:
        f.write(rendered_train_config)

    with open(embed_config_path, "w") as f:
        f.write(rendered_embed_config)

    with open(register_config_path, "w") as f:
        f.write(rendered_register_config)

    return train_config_path, embed_config_path, register_config_path


def create_pca_configs_from_yaml(config_yaml_path: str) -> t.List[t.Dict[str, str]]:
    with open(config_yaml_path, "r") as file:
        data = yaml.safe_load(file)

    pipeline_configs = data["pipeline_configs"]
    config_paths = []

    for pipeline_config in pipeline_configs:
        _validate_train_embed_config(pipeline_config)
        train_config_path, embed_config_path, register_model_config_path = create_train_embed_configs(**pipeline_config)
        config_paths.append(
            {
                "train_gcs_config_path": train_config_path,
                "embed_gcs_config_path": embed_config_path,
                "register_model_gcs_config_path": register_model_config_path,
            }
        )

    return config_paths
