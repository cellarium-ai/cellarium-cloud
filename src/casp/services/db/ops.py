import typing as t

from casp.services.db import db_session

if t.TYPE_CHECKING:
    from casp.services.db import models


def increment_user_cells_processed(user: models.User, number_of_cells: int) -> None:
    """
    Increment `requests_processed` and `cells_processed user stats`.
    :param user: A SQLAlchemy user model
    :param number_of_cells: number of cells to increment by
    """
    user.requests_processed += 1
    user.cells_processed += number_of_cells
    db_session.commit()