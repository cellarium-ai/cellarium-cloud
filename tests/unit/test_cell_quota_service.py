"""
Tests some methods in the cellarium_general data manager
"""

import datetime

from google.cloud import bigquery
from mockito import unstub, when

from casp.services import db, utils
from casp.services.api import schemas
from casp.services.api.services.cell_quota_service import CellQuotaService
from casp.services.db.models.users import User


class TestCellQuotaService:

    def setup_method(self) -> None:
        # Gotta mock up the google service stuff that gets called in the constructor so the tests
        # don't fail when you run them on a system without credentials set up
        # Important note: Any tests that actually rely on bigquery or the db will fail, so the
        # tests I've written here mock up calls to methods that need those
        when(utils).get_google_service_credentials().thenReturn((None, None))
        when(bigquery).Client(credentials=None, project=None).thenReturn(None)
        when(db).get_db_session_maker().thenReturn(lambda: None)
        self.cell_quota_service = CellQuotaService()

    def teardown_method(self) -> None:
        unstub()

    def test_get_quota_for_user(self):
        user = User(
            id=1,
            cell_quota=100000,
            lifetime_cell_quota=200000,
            quota_increased=True,
            total_cells_processed=75000,
            created_at=datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0),
        )
        when(self.cell_quota_service.cell_quota_dm).get_cells_processed_this_week_for_user(user=user).thenReturn(75000)

        assert self.cell_quota_service.get_quota_for_user(user=user) == schemas.UserQuota(
            user_id=1,
            weekly_quota=100000,
            remaining_weekly_quota=25000,
            quota_reset_date=user.created_at + datetime.timedelta(days=7),
            lifetime_quota=200000,
            remaining_lifetime_quota=125000,
            quota_increased=True
        )

    def test_get_quota_for_user_lifetime_low(self):
        user = User(
            id=1,
            cell_quota=100000,
            lifetime_cell_quota=200000,
            quota_increased=True,
            total_cells_processed=190000,
            created_at=datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0),
        )
        when(self.cell_quota_service.cell_quota_dm).get_cells_processed_this_week_for_user(user=user).thenReturn(0)

        assert self.cell_quota_service.get_quota_for_user(user=user) == schemas.UserQuota(
            user_id=1,
            weekly_quota=100000,
            remaining_weekly_quota=10000,
            quota_reset_date=user.created_at + datetime.timedelta(days=7),
            lifetime_quota=200000,
            remaining_lifetime_quota=10000,
            quota_increased=True
        )
