from google.cloud import bigquery

from casp.data_manager import cache
from casp.services import utils
from casp.services.db import init_db


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
    redis_cache_key_prefix: str

    def __init__(self):
        credentials, project = utils.get_google_service_credentials()
        self.credentials = credentials
        self.project = project
        self.block_coo_matrix_db_client = bigquery.Client(credentials=self.credentials, project=self.project)
        self.system_data_db_session = init_db()
        self.cache = cache.RedisCache(cache_key_prefix=self.redis_cache_key_prefix)
