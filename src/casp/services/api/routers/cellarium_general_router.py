import typing as t

from fastapi import APIRouter, Depends, Security

from casp.services.api import dependencies, schemas, services
from casp.services.db import models

cellarium_general_router = APIRouter(prefix="/cellarium-general")
cellarium_general_service = services.CellariumGeneralService()


@cellarium_general_router.get("/application-info", response_model=schemas.ApplicationInfo)
async def application_info():
    """
    Get Cellarium CAS application info such as version, default feature schema name.
    """
    return cellarium_general_service.get_application_info()


@cellarium_general_router.get("/validate-token")
async def validate_token(_: models.User = Depends(dependencies.authenticate_user)):
    """
    Validate authorization token from `Bearer` header

    :return: Success message if token is valid, otherwise return 401 Unauthorized status code if token
    is invalid or missing
    """
    return {"detail": "Success"}


@cellarium_general_router.get("/feature-schemas", response_model=t.List[schemas.FeatureSchemaInfo])
async def get_feature_schemas(_: models.User = Depends(dependencies.authenticate_user)):
    """
    Get list of all Cellarium CAS feature schemas

    :return: List of feature schema info objects
    """
    return cellarium_general_service.get_feature_schemas()


@cellarium_general_router.get("/feature-schema/{schema_name}", response_model=t.List[str])
async def get_feature_schema_by(schema_name: str, _: models.User = Depends(dependencies.authenticate_user)):
    """
    Get a specific feature schema by its unique name

    :param schema_name: unique feature schema name

    :return: List of features in a correct order
    """
    return cellarium_general_service.get_feature_schema_by(schema_name=schema_name)


@cellarium_general_router.get("/list-models", response_model=t.List[schemas.CASModel])
async def get_model_list(user: models.User = Depends(dependencies.authenticate_user)):
    """
    Get list of all Cellarium CAS models for a specific user based on their permissions.

    :return: List of CAS models
    """
    return cellarium_general_service.get_model_list_for_user(user=user)
