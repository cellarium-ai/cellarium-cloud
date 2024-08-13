from kfp import dsl


# @dsl.component()
def mean_var_std_train(gcs_config_path: str):
    """
    Cellarium ML CLI component for running one pass mean var std training.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK") is not None:
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["onepass_mean_var_std", "fit", "--config", gcs_config_path])


# @dsl.component()
def tdigest_train(gcs_config_path: str):
    """
    Cellarium ML CLI component for running tdigest model training.

    :param gcs_config_path: Config file path on GCS.
    """
    import os

    from cellarium.ml.cli import main as cellarium_ml_cli

    if os.environ.get("RANK") is not None:
        os.environ["NODE_RANK"] = os.environ.get("RANK")

    cellarium_ml_cli(args=["tdigest", "fit", "--config", gcs_config_path])


# @dsl.component()
def tdigest_filter_features_component(gcs_config_path: str):
    """
    Get TDigest model and expedite its `median_g` argument. Create a new feature set with only the features that
    don't have a median_g value of nan.

    :param gcs_config_path: Config file path on GCS.
    """
    import io

    import pandas as pd
    import yaml
    from cellarium.ml import CellariumModule
    from smart_open import open

    with open(gcs_config_path, "r") as file:
        config_data = yaml.safe_load(file)

    tdigest_model_path = config_data["tdigest_model_path"]
    feature_filter_path = config_data["feature_file_filter_path"]

    with open(tdigest_model_path, "rb") as f:
        ckpt = io.BytesIO(f.read())

    module = CellariumModule.load_from_checkpoint(ckpt)
    median_g_non_mask = ~module.model.median_g.isnan()
    new_var_names = module.model.var_names_g[median_g_non_mask]

    df = pd.DataFrame({"original_feature_id": new_var_names})
    df.to_csv(feature_filter_path)
