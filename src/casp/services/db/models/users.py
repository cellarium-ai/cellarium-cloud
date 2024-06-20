import datetime
import enum

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import backref, column_property, deferred, relationship

from casp.services import db, settings

# The current version of the token
CURRENT_TOKEN_VERSION: str = "2"


class User(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    email = sa.Column(sa.String(255), unique=True, nullable=False)
    username = sa.Column(sa.String(255), nullable=False)
    active = sa.Column(sa.Boolean(), default=True, nullable=False)
    requests_processed = sa.Column(sa.Integer, default=0, nullable=False)
    cells_processed = sa.Column(sa.Integer, default=0, nullable=False)
    is_admin = sa.Column(sa.Boolean(), default=True, nullable=False)
    cell_quota = sa.Column(sa.Integer(), default=50000, nullable=False)
    created_at = sa.Column(sa.DateTime, nullable=False, default=datetime.datetime.now(datetime.timezone.utc))
    ask_for_feedback = sa.Column(sa.Boolean(), default=True, nullable=False)

    __tablename__ = "users_user"

    def __repr__(self):
        return self.email


class UserActivityEvent(enum.Enum):
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    STARTED = "STARTED"


class UserActivity(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    user_id = sa.Column(sa.Integer, sa.ForeignKey(f"{User.__tablename__}.id"))
    request_id = sa.Column(sa.String(255), nullable=True)
    cell_count = sa.Column(sa.Integer, default=0, nullable=False)
    model_name = sa.Column(sa.String(255), nullable=False)
    method = sa.Column(sa.String(255), nullable=True)
    finished_time = sa.Column(sa.DateTime, default=datetime.datetime.now(datetime.timezone.utc))
    event = sa.Column(sa.Enum(UserActivityEvent), nullable=False)

    __tablename__ = "users_useractivity"


# Add properties to user model for metrics
sa.inspect(User).add_property(
    key="total_cells_processed",
    prop=column_property(
        User.cells_processed
        + sa.select(sa.func.sum(UserActivity.cell_count))
        .where((UserActivity.user_id == User.id) & (UserActivity.event == UserActivityEvent.SUCCEEDED))
        .scalar_subquery()
    ),
)
sa.inspect(User).add_property(
    key="total_requests_processed",
    prop=column_property(
        User.requests_processed
        + sa.select(sa.func.count(UserActivity.id))
        .where((UserActivity.user_id == User.id) & (UserActivity.event == UserActivityEvent.SUCCEEDED))
        .scalar_subquery()
    ),
)


class UserKey(db.Base):
    id = sa.Column(sa.Integer, primary_key=True)
    key_locator = sa.Column(UUID, nullable=False, unique=True)
    user_id = sa.Column(sa.Integer, sa.ForeignKey(f"{User.__tablename__}.id"), nullable=False)
    created_date = sa.Column(sa.DateTime, default=datetime.datetime.now(datetime.timezone.utc), nullable=False)
    active = sa.Column(sa.Boolean(), default=True, nullable=False)
    expires = sa.Column(
        sa.DateTime,
        default=datetime.datetime.now(datetime.timezone.utc)
        + datetime.timedelta(seconds=settings.JWT_DEFAULT_TOKEN_TTL),
        nullable=False,
    )

    # A hash of the user's key.  Adding deferred to prevent loading the key unless it is needed to avoid leaking it accidentally
    key_hash = deferred(sa.Column(sa.Text, nullable=False))

    # Causes the user to be fetched when the key is fetched which makes verification much easier
    user = relationship("User", backref=backref("user_keys", uselist=True, cascade="all, delete-orphan"))

    __tablename__ = "users_userkey"
