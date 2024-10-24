import typing as t

from fastapi import APIRouter, Depends
from fastapi.responses import RedirectResponse

from casp.services import constants, settings
from casp.services.api import dependencies, schemas, services
from casp.services.db import models

cellarium_general_router = APIRouter(prefix="/cellarium-general")
cellarium_general_service = services.CellariumGeneralService()
cell_quota_service = services.CellQuotaService()


@cellarium_general_router.get("/application-info", response_model=schemas.ApplicationInfo)
async def application_info():
    """
    Get Cellarium CAS application info such as version, default feature schema name.
    """
    return cellarium_general_service.get_application_info()


@cellarium_general_router.post("/validate-client-version", response_model=schemas.ClientVersionOutput)
async def validate_client_version(client_version_info: schemas.ClientVersionInput):
    """
    Check whether the client version is new enough to work with the server
    """
    return schemas.ClientVersionOutput(
        is_valid=cellarium_general_service.validate_client_version(client_version=client_version_info.client_version),
        min_version=constants.MIN_CLIENT_VERSION,
    )


@cellarium_general_router.get("/validate-token", response_model=schemas.UserInfo)
async def validate_token(user: models.User = Depends(dependencies.authenticate_user)):
    """
    Validate authorization token from `Bearer` header

    :return: User information if the token is valid, otherwise return 401 Unauthorized status code if token
    is invalid or missing
    """
    return schemas.UserInfo(username=user.username, email=user.email, should_ask_for_feedback=user.ask_for_feedback)


@cellarium_general_router.get("/feedback/answer")
async def feedback_answer(client_session_id: str, client_action_id: str):
    """
    Redirect to the feedback form with the client session ID
    """
    return RedirectResponse(
        settings.FEEDBACK_FORM_BASE_URL.format(client_session_id=client_session_id, client_action_id=client_action_id)
    )


@cellarium_general_router.post("/feedback/opt-out", response_model=schemas.UserInfo)
async def feedback_opt_out(user: models.User = Depends(dependencies.authenticate_user)):
    """
    Allow user to opt out of feedback requests

    :return: User information if reflecting the updated feedback preference
    """
    updated_user = cellarium_general_service.feedback_opt_out(user=user)
    return schemas.UserInfo(
        username=updated_user.username, email=updated_user.email, should_ask_for_feedback=updated_user.ask_for_feedback
    )


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


@cellarium_general_router.get("/quota", response_model=schemas.UserQuota)
async def get_user_quota(user: models.User = Depends(dependencies.authenticate_user)):
    """
    Get user quota information

    :return: User quota information
    """
    return cell_quota_service.get_quota_for_user(user=user)

@cellarium_general_router.post("/increase-quota/{user_email}")
async def increase_quota_for_user_by_email(
    user_email: str,
    admin_user: models.User = Depends(dependencies.authenticate_user)
):
    """
    Increase the lifetime cell quota for the user specified by user_email if their lifetime quota
    has not yet been increased

    :param user_email: Email of the user whose quota will be increased
    """
    user_for_increase = cellarium_general_service.get_user_by_email(user_email=user_email)
    cell_quota_service.increase_quota(admin_user=admin_user, user_for_increase=user_for_increase)
