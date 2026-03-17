from cellarium.cas_backend.core import config

if config.ENV_TYPE == "local":
    settings = config.LocalSettings()
elif config.ENV_TYPE == "production":
    settings = config.ProductionSettings()
elif config.ENV_TYPE == "test":
    settings = config.TestSettings()
else:
    settings = config.DevSettings()
