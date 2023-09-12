import typing as t

from casp.services.db import db_session, models


def increment_user_cells_processed(user: models.User, number_of_cells: int) -> None:
    """
    Increment `requests_processed` and `cells_processed user stats`.
    :param user: A SQLAlchemy user model
    :param number_of_cells: number of cells to increment by
    """
    user.requests_processed += 1
    user.cells_processed += number_of_cells
    db_session.commit()


def get_models_for_user(user: "models.User") -> t.List[models.CASModel]:
    """
    Retrieve query with available models for a specific user depending on their status
    Admin users would get all models while non-admin status users would get admin_use_only
    :param user: CAS user model instance
    :return: Queryset with CAS models
    """
    if not user.is_admin:
        # Return only selected models for non-admin users
        return models.CASModel.query.filter_by(admin_use_only=False).all()

    return models.CASModel.query.all()


def get_model_by(model_name: str) -> models.CASModel:
    """
    Retrieve CAS model by its system name
    :param model_name: Model system name
    :return: CAS model
    """
    return models.CASModel.query.filter_by(model_name=model_name).first()


def set_default_model_by(model_id: int) -> None:
    """
    Change system's default CAS model. Changes is_default_model by id and sets all other models being default to False
    :param model_id: ID of the model
    """
    models.CASModel.query.update({models.CASModel.is_default_model: False})
    db_session.commit()
    models.CASModel.query.filter_by(id=model_id).update({models.CASModel.is_default_model: True})
    db_session.commit()
