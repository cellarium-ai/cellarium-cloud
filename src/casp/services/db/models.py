import sqlalchemy as sa

from casp.services import db


class User(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    email = sa.Column(sa.String(255), unique=True, nullable=False)
    username = sa.Column(sa.String(255), nullable=False)
    active = sa.Column(sa.Boolean(), default=True, nullable=False)
    cas_request_count = sa.Column(sa.Integer, default=0, nullable=False)
    cas_scRNA_cells_processed = sa.Column(sa.Integer, default=0, nullable=False)

    __tablename__ = "user"

    def __repr__(self):
        return self.email
