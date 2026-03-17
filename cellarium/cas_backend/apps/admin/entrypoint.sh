#!/bin/bash
# Update Database Migrations
alembic -c casp/services/db/alembic.ini upgrade head
# Run Flask App with Gunicorn
gunicorn --bind 0.0.0.0:8000 --workers 1 --threads 8 --timeout 0 casp.services.admin.server:flask_app