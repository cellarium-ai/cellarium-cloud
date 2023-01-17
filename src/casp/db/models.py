import sqlalchemy as sa
from casp import db


class User(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    email = sa.Column(sa.String(255), unique=True)
    username = sa.Column(sa.String(255))
    active = sa.Column(sa.Boolean())
    cas_request_count = sa.Column(sa.Integer)
    cas_scRNA_cells_processed = sa.Column(sa.Integer)

    __tablename__ = 'user'

    def __repr__(self):
        return self.email
