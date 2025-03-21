import typing as t

import sqlalchemy as sa
from google.cloud import bigquery

from casp.services import utils
from casp.services.db import get_db_session_maker


class BaseDataManager:
    """
    Base class for all Data Managers.  It also provides a BigQuery client and a Postgres DB session.
    Data Managers are used by services to communicate with databases.

    Class attributes:
    credentials: Google service credentials object
    project: Google Cloud project ID
    bigquery_client: Google BigQuery client object to communicate with BigQuery
    postgres_db_session: Postgres DB session object to communicate with Postgres DB
    SQL_TEMPLATE_DIR (has to be set): Path to a directory with SQL templates related to Cellarium API.
    """

    SQL_TEMPLATE_DIR: str

    def __init__(self):
        credentials, project = utils.get_google_service_credentials()
        self.credentials = credentials
        self.project = project
        self.block_coo_matrix_db_client = bigquery.Client(credentials=self.credentials, project=self.project)
        self.system_data_db_session_maker = get_db_session_maker()

    @staticmethod
    def batch_bulk_insert(
        table: sa.Table,
        rows_to_insert: t.List[t.Dict[str, t.Any]],
        connection: sa.engine.Connection,
        batch_size: int = 10000,
    ):
        for i in range(0, len(rows_to_insert), batch_size):
            batch = rows_to_insert[i : i + batch_size]
            connection.execute(table.insert().values(batch))
