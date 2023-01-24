from casp import _settings

if _settings.ENV_TYPE == "local":
    settings = _settings.LocalSettings()
elif _settings.ENV_TYPE == "development":
    settings = _settings.DevSettings()
elif _settings.ENV_TYPE == "production":
    settings = _settings.ProductionSettings()
else:
    raise ValueError(
        f"Couldn't parse environment type. It should be one of: `local`, `development`, `production`"
        f"Instead it got: `{_settings.ENV_TYPE}`"
    )
