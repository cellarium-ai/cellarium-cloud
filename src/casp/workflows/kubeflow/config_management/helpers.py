import json
import typing as t

import yaml
from mako.template import Template
from smart_open import open

from casp.services import settings

PIPELINE_CONFIGS_TEMPLATE_PATH = f"{settings.CAS_DIR}/workflows/kubeflow/config_management/config_templates"
CONFIG_BUCKET_PATH = "gs://cellarium-file-system/ml-configs/auto-generated"


def get_component_config_template_kwargs(config_data: t.Dict[str, t.Any], component_name: str) -> t.Dict[str, t.Any]:
    """
    Extract the config data for a specific component config from the config data dictionary.

    It will produce a dictionary that can be used for template rendering. To define the keys that are specific to a
    component, the key in the config data dictionary should be in the format of
    <component_name>__<component_specific_key>.
    It first looks at the keys that are shared, then looks at the keys that are specific to the component. If there are
    any shared keys that are also specific to the component, the component specific key will overwrite the shared key.
    Note: the result dictionary doesn't have information about the config template variables, so in case if a shared key
    is not present in the template, it still will be included in `result` dictionary. However, template generator is
    smart enough to ignore the extra keys.

    :param config_data: Dictionary of config data for a specific pipeline.
    :param component_name: Name of the component to extract config data for.

    :return: Dictionary of config data for a specific component's config template.
    """
    component_related_keys = [key for key in config_data.keys() if component_name in key]
    shared_keys = [key for key in config_data.keys() if "__" not in key]

    result = {key: config_data[key] for key in shared_keys}

    for key in component_related_keys:
        result[key.split("__")[1]] = config_data[key]

    return result


def get_component_names_from_config_data(config_data: t.Dict[str, t.Any]) -> t.List[str]:
    """
    Extract the component names from the config data dictionary. Look at the config data keys and extract the component
    names from the keys that are in the format of <component_name>__<component_specific_key>.

    :param config_data: Dictionary of config data for a specific pipeline.

    :return: List of component names.
    """
    component_names = set()

    for key in config_data.keys():
        if "__" in key:
            component_name = key.split("__")[0]
            component_names.add(component_name)

    return list(component_names)


def render_and_save_template(template_path: str, context: dict, output_path: str):
    """
    Render a Mako template with the provided context and save the output to a file.

    :param template_path: Path to the Mako template file.
    :param context: Dictionary containing the context to render the template.
    :param output_path: Path where the rendered output should be saved (can be either local or GCS path).
    """
    template = Template(filename=template_path)
    rendered_content = template.render(**context)

    with open(output_path, "w") as f:
        f.write(rendered_content)


def create_configs(config_yaml_path: str) -> t.List[str]:
    """
    Create and save the pipeline config files based on the config data in the provided YAML file.

    YAML file should have a key called `pipeline_configs` which is a list of objects. Each object should consist of
    shared config data and component specific config data. Shared config keys are all keys that don't have a component
    name in them. Component specific config keys are all keys that have a component name in them following the format:
    <component_name>__<component_specific_key>. For example, `pca_train__config_unique_name`. There is one required
    component specific key for each component: <component_name>__config_unique_name. This key is used to generate the
    config file name. For example, `pca_train__config_unique_name: pca_train_1` will generate a config file with the
    name of `pca_train-pca_train_1_config.yml`. Missing this key will result in a missing config file.

    :param config_yaml_path: Path to the YAML file containing the config data.
    """
    with open(config_yaml_path, "r") as file:
        data = yaml.safe_load(file)

    pipeline_configs = data["pipeline_configs"]

    print("Creating config files")
    config_paths = []
    for pipeline_config in pipeline_configs:
        curriculum_name = pipeline_config["curriculum_name"]

        component_names = get_component_names_from_config_data(pipeline_config)
        config_object = {}

        for component_name in component_names:
            template_path = f"{PIPELINE_CONFIGS_TEMPLATE_PATH}/{component_name}_config.yml.mako"
            component_config_unique_name = pipeline_config[f"{component_name}__config_unique_name"]
            component_config_filename = f"{component_name}-{component_config_unique_name}_config.yml"
            output_path = f"{CONFIG_BUCKET_PATH}/{curriculum_name}/{component_config_filename}"

            context = get_component_config_template_kwargs(config_data=pipeline_config, component_name=component_name)
            render_and_save_template(template_path=template_path, context=context, output_path=output_path)

            config_object[f"{component_name}_gcs_config_path"] = output_path

        config_paths.append(config_object)
    print(f"Config files created and saved in {CONFIG_BUCKET_PATH}")

    # DSL Pipelines require JSON strings as input
    config_paths = [json.dumps(config) for config in config_paths]

    return config_paths