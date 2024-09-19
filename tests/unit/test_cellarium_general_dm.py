"""
Tests some methods in the cellarium_general data manager
"""

import pytest
from mockito import unstub, when

from casp.services.api.data_manager.cellarium_general import CellariumGeneralDataManager
from casp.services.db.models.users import User

class TestCellariumGeneralDM:

    def setup_method(self) -> None:
        self.cellarium_general_dm = CellariumGeneralDataManager()

    def teardown_method(self) -> None:
        unstub()

    def test_get_remaining_quota_for_user_some_left(self):
        user = User(
            id = 1,
            cell_quota = 100000,
            lifetime_cell_quota = 200000,
            total_cells_processed = 0
        )
        when(self.cellarium_general_dm).get_cells_processed_this_week_for_user(user).thenReturn(25000)
        
        assert self.cellarium_general_dm.get_remaining_quota_for_user(user=user) == 75000

    def test_get_remaining_quota_for_user_none_left_weekly(self):
        user = User(
            id = 1,
            cell_quota = 100000,
            lifetime_cell_quota = None,
            total_cells_processed = 0
        )
        when(self.cellarium_general_dm).get_cells_processed_this_week_for_user(user).thenReturn(100001)
        
        assert self.cellarium_general_dm.get_remaining_quota_for_user(user=user) == 0

    def test_get_remaining_quota_for_user_none_left_lifetime(self):
        user = User(
            id = 1,
            cell_quota = 100000,
            lifetime_cell_quota = 200000,
            total_cells_processed = 200000
        )
        when(self.cellarium_general_dm).get_cells_processed_this_week_for_user(user).thenReturn(0)
        
        assert self.cellarium_general_dm.get_remaining_quota_for_user(user=user) == 0
    