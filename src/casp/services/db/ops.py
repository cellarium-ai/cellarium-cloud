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
    if not user.is_admin:
        # Return only selected models for non-admin users
        return models.CASModel.query.filter_by(admin_use_only=False).all()

    return models.CASModel.query.all()


def get_model_by(system_name: str) -> models.CASModel:
    return models.CASModel.query.filter_by(system_name=system_name).first()
