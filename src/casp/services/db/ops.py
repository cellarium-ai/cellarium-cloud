import typing as t
from casp.services.db import db_session

if t.TYPE_CHECKING:
    from casp.services.db import models


def increment_user_cells_processed(user: "models.User", number_of_cells: int):
    user.cas_request_count += 1
    user.cas_scRNA_cells_processed += number_of_cells
    db_session.commit()
