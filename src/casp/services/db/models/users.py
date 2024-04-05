import datetime

import sqlalchemy as sa
from sqlalchemy.orm import column_property

from casp.services import db


class User(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    email = sa.Column(sa.String(255), unique=True, nullable=False)
    username = sa.Column(sa.String(255), nullable=False)
    active = sa.Column(sa.Boolean(), default=True, nullable=False)
    requests_processed = sa.Column(sa.Integer, default=0, nullable=False)
    cells_processed = sa.Column(sa.Integer, default=0, nullable=False)
    is_admin = sa.Column(sa.Boolean(), default=True, nullable=False)

    __tablename__ = "users_user"

    def __repr__(self):
        return self.email


class UserActivity(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    user_id = sa.Column(sa.Integer, sa.ForeignKey(f"{User.__tablename__}.id"))
    cell_count = sa.Column(sa.Integer, default=0, nullable=False)
    model_name = sa.Column(sa.String(255), nullable=False)
    method = sa.Column(sa.String(255), nullable=True)
    finished_time = sa.Column(sa.DateTime, default=datetime.datetime.now(datetime.timezone.utc))

    __tablename__ = "users_useractivity"


# Add properties to user model for metrics
sa.inspect(User).add_property(
    key="total_cells_processed",
    prop=column_property(
        User.cells_processed
        + sa.select(sa.func.sum(UserActivity.cell_count)).where(UserActivity.user_id == User.id).scalar_subquery()
    ),
)
sa.inspect(User).add_property(
    key="total_requests_processed",
    prop=column_property(
        User.requests_processed
        + sa.select(sa.func.count(UserActivity.id)).where(UserActivity.user_id == User.id).scalar_subquery()
    ),
)
