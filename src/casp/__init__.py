from casp import _settings

if _settings.ENV_TYPE == "local":
    settings = _settings.LocalSettings()
elif _settings.ENV_TYPE == "production":
    settings = _settings.ProductionSettings()
else:
    settings = _settings.DevSettings()
