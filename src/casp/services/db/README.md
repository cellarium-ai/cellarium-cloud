# CAS Database Service
This module requires running database cluster.
Create db cluster:
```
gcloud sql instances create cas-db \
--database-version=POSTGRES_14 \
--cpu=2 \
--memory=7680MB \
--network=ai-matching \
--region=us-central
```
Create a cas-db-user:
```
gcloud sql users create cas-db-user \
--instance=cas-db \
--password={password-goes-here}
```
Create a cas-db:
```
gcloud sql databases create cas-db \
--instance=cas-db \

```
## Code Base Info

Database module `db`consist of:
- `migrations` Database migration history goes here;
- `models.py` All database models go here (please, feel free to extend it to module if necessary);
- `alembic.ini` Provides metadata for migration manager;

## Database Adapter and ORM
[SQLAlchemy](https://www.sqlalchemy.org/) is used for db object-relational mapping. \
[Psycopg](https://pypi.org/project/psycopg2/) is used as a default database adapter. \
`DB_HOST`, `DB_PORT`, `DB_USER`, `DB_PASSWORD`, `DB_NAME` are required secret environment variables that are used to configure `SQLALCHEMY_DATABASE_URI` in `casp/settings.py`.
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
