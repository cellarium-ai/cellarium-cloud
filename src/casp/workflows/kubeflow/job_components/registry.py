from kfp import dsl


@dsl.component()
def register_embedding_model(gcs_config_path: str) -> None:
    """
    Component for registering embedding model in the registry by creating an instance in the database with the model
    information

    :param gcs_config_path: GCS path to the config file with the model information
    """
    import yaml
    from smart_open import open

    from casp.workflows.kubeflow.job_components_library.clients import RegistryClient

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    client = RegistryClient()
    client.register_model(
        model_name=config_data["model_name"],
        model_file_path=config_data["model_file_path"],
        embedding_dimension=config_data["embedding_dimension"],
        bq_dataset_name=config_data["bq_dataset_name"],
        schema_name=config_data["schema_name"],
    )
