import datetime

from casp.services.api import schemas
from casp.services.api.data_manager import CellQuotaDataManager
from casp.services.db import models


class CellQuotaService:
    """
    Cell Quota Service. Service providing functions for getting info on user cell quotas.
    """

    def __init__(self, cell_quota_dm: CellQuotaDataManager = None):
        self.cell_quota_dm = cell_quota_dm or CellQuotaDataManager()

    def get_remaining_quota_for_user(self, user: models.User) -> int:
        """
        Get remaining cell processing quota for a specific user (i.e. their quota minus cells processed this week
        or their remaining lifetime quota if they have one)

        :param user: User object to check quota for

        :return: Remaining cell processing quota
        """
        # Return the smaller of the remaining weekly and lifetime quotas, or 0 if it's negative
        # (which might happen if the user exceeds their quota, either because we allow them to, or
        # because of any weirdness in the simultaneous processing of cells)
        weekly_remaining_quota = user.cell_quota - self.cell_quota_dm.get_cells_processed_this_week_for_user(user=user)

        lifetime_remaining_quota = (
            user.lifetime_cell_quota - user.total_cells_processed
            if user.lifetime_cell_quota is not None
            else float("inf")
        )

        remaining_quota = min(weekly_remaining_quota, lifetime_remaining_quota)

        # Ensure the remaining quota is not negative
        return max(remaining_quota, 0)

    def get_quota_for_user(self, user: models.User) -> schemas.UserQuota:
        """
        Get user quota information

        :param user: User object to check quota for

        :return: User quota information
        """

        days_to_add = datetime.datetime.now().weekday() - user.created_at.weekday()
        if days_to_add <= 0:
            days_to_add += 7
        quota_reset_date = datetime.datetime.today().replace(
            hour=0, minute=0, second=0, microsecond=0
        ) + datetime.timedelta(days=days_to_add)

        remaining_lifetime_quota = (
            user.lifetime_cell_quota - user.total_cells_processed if user.lifetime_cell_quota is not None else None
        )

        return schemas.UserQuota(
            user_id=user.id,
            weekly_quota=user.cell_quota,
            remaining_weekly_quota=self.get_remaining_quota_for_user(user=user),
            quota_reset_date=quota_reset_date,
            lifetime_quota=user.lifetime_cell_quota,
            remaining_lifetime_quota=remaining_lifetime_quota,
        )
