# CAS Database Service
This module requires running database cluster.
Create a db cluster:
```
gcloud sql instances create cas-db-cluster \
--database-version=POSTGRES_14 \
--cpu=1 \
--memory 3.75GB \
--storage-size 10GB \
--require-ssl \
--region=us-central
```
Create a cas-db-user:
```
gcloud sql users create cas-db-user \
--instance=cas-db-cluster \
--password=b407f042a92dfac5913764ab
```
Create a cas-db:
```
gcloud sql databases create cas-db \
--instance=cas-db-cluster
```
## Code Base Info

Database module `db`consist of:
- `migrations` Database migration history goes here;
- `models.py` All database models go here;
- `ops.py` Data operations go here;
- `alembic.ini` Provides metadata for migration manager;

## Database Adapter and ORM
[SQLAlchemy](https://www.sqlalchemy.org/) is used for db object-relational mapping. \
[pg8000](https://pypi.org/project/pg8000/) is used as a database adapter. 
## Environment Variables and Connection
Locally connection is approached through regular postgresql connection (by port). In the cloud it's done through proxy. To connect to proxy it uses a unix socket. \
`SQLALCHEMY_DATABASE_URI`  is used by [SQLAlchemy](https://www.sqlalchemy.org/) to connect and configured from the following secret environment variables:
* `DB_HOST`, `DB_PORT`, `DB_USER`, `DB_PASSWORD`, `DB_NAME` in local environment.
* `DB_USER`, `DB_PASSWORD`, `DB_NAME`, `DB_INSTANCE_UNIX_SOCKET` in development and production environments. 
## Database Migrations
[Alembic](https://alembic.sqlalchemy.org/en/latest/) is used for managing database migrations. \
Each time CAS Database models are updated it is required to:
Generate a new migration: 
```
alembic -c casp/services/db/alembic.ini revision --autogenerate -m "{migration-message-goes-here}"
```
Apply migrations to the database:
```
alembic -c casp/services/db/alembic.ini upgrade head 
```
