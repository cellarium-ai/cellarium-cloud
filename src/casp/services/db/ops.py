from casp.services.db import db_session, models


def set_default_model_by(model_id: int) -> None:
    """
    Change system's default CAS model. Changes is_default_model by id and sets all other models being default to False
    :param model_id: ID of the model
    """
    models.CASModel.query.update({models.CASModel.is_default_model: False})
    db_session.commit()
    models.CASModel.query.filter_by(id=model_id).update({models.CASModel.is_default_model: True})
    db_session.commit()
