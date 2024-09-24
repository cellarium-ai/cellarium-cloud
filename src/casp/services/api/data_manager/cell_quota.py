import datetime

import sqlalchemy as sa

from casp.data_manager import BaseDataManager
from casp.services import settings
from casp.services.db import models


class CellQuotaDataManager(BaseDataManager):

    # SQL query directory
    SQL_QUERY_DIR = f"{settings.SERVICES_DIR}/api/data_manager/sql_queries"

    # SQL queries from files
    with open(f"{SQL_QUERY_DIR}/get_cells_processed_this_week_for_user.sql.sa", "r") as f:
        SQL_GET_CELLS_PROCESSED_THIS_WEEK = f.read()

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
        weekly_remaining_quota = user.cell_quota - self.get_cells_processed_this_week_for_user(user=user)

        lifetime_remaining_quota = (
            user.lifetime_cell_quota - user.total_cells_processed
            if user.lifetime_cell_quota is not None
            else float("inf")
        )

        remaining_quota = min(weekly_remaining_quota, lifetime_remaining_quota)

        # Ensure the remaining quota is not negative
        return max(remaining_quota, 0)

    def get_cells_processed_this_week_for_user(self, user: models.User) -> int:
        """
        Get number of cells processed by a specific user this week

        :param user: User object to check cells processed for

        :return: Number of cells processed this week
        """
        # Start of week for a user is based on the day of the week they signed up, so we need to
        # use the day of the week they signed up to calculate the start of the week
        current_weekday = datetime.datetime.today().weekday()
        start_of_the_week_for_user = user.created_at.weekday()
        days_to_subtract = current_weekday - start_of_the_week_for_user
        if days_to_subtract < 0:
            days_to_subtract += 7
        start_of_week = datetime.datetime.today().replace(
            hour=0, minute=0, second=0, microsecond=0
        ) - datetime.timedelta(days=days_to_subtract)

        with self.system_data_db_session_maker() as session:
            cell_count_sum_query = sa.sql.text(self.SQL_GET_CELLS_PROCESSED_THIS_WEEK)
            return session.execute(cell_count_sum_query, {"id": user.id, "start_of_week": start_of_week}).first()[0]
