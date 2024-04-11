from casp.services.db import get_db_session_maker, models


def set_default_model_by(model_id: int) -> None:
    """
    Change system's default CAS model. Changes is_default_model by id and sets all other models being default to False
    :param model_id: ID of the model
    """
    with get_db_session_maker()() as db_session:
        with db_session.begin():
            db_session.query(models.CASModel).update({models.CASModel.is_default_model: False})
            db_session.query(models.CASModel).filter_by(id=model_id).update({models.CASModel.is_default_model: True})
