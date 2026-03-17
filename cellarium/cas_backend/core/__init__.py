from cellarium.cas_backend.core import _settings

if _settings.ENV_TYPE == "local":
    settings = _settings.LocalSettings()
elif _settings.ENV_TYPE == "production":
    settings = _settings.ProductionSettings()
elif _settings.ENV_TYPE == "test":
    settings = _settings.TestSettings()
else:
    settings = _settings.DevSettings()
