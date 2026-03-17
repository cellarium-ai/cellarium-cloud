# Project Refactoring Summary

**Date**: 2026-03-17
**Branch**: `fg/project-structure-refactor`

## ✅ Completed

### 1. Package Restructuring
- ✅ Renamed `casp` → `cellarium.cas_backend` (namespace package)
- ✅ Removed `src/` wrapper
- ✅ Created clean architecture:
  - `cellarium/cas_backend/apps/` - Application services
  - `cellarium/cas_backend/core/` - Shared infrastructure
- ✅ Renamed `api` → `compute` app
- ✅ Unified data_managers (merged from 2 locations)
- ✅ Moved `db/` to core (not app-level)
- ✅ Removed underscore prefixes (`_auth` → `auth`, `_settings.py` → `config.py`)

### 2. Deployment Restructuring
- ✅ Created `deploy/` directory at repo root
- ✅ Moved all Dockerfiles to `deploy/docker/`
- ✅ Moved deployment scripts to `deploy/scripts/`
- ✅ Moved CloudRun configs to `deploy/cloudrun/`

### 3. Scripts & Workflows
- ✅ Moved `scripts/` to repo root
- ✅ Moved `workflows/` to repo root (kubeflow, wdl)
- ✅ Created `scripts/bq_ops/entrypoints/` for WDL entry points

### 4. Dependency Management
- ✅ Created `pyproject.toml` with Poetry
- ✅ Generated `poetry.lock` for reproducible builds
- ✅ Deleted old requirements files (requirements.txt, dev-requirements.txt, etc.)
- ✅ Auto-generate requirements.txt from Poetry via `make requirements`
- ✅ Added Makefile for common tasks

### 5. Docker Updates
- ✅ Created new Dockerfiles using new structure:
  - `Dockerfile.compute` - Main API service
  - `Dockerfile.admin` - Admin interface
  - `Dockerfile.model_inference` - ML serving
- ✅ Updated to copy `cellarium/` instead of `src/`
- ✅ Updated PYTHONPATH for namespace package
- ✅ Comprehensive `.dockerignore`
- ✅ Docker build documentation

### 6. Import Updates
- ✅ Auto-updated 92 Python files with new import paths
- ✅ Created `update_imports.py` migration script
- ✅ All imports: `casp.*` → `cellarium.cas_backend.*`

## 📋 TODO

### Critical (Blocking)

1. **Update tox.ini**
   - Change `--cov=casp` → `--cov=cellarium.cas_backend`
   - Update paths in lint/format commands

2. **Update WDL workflows**
   - Find/replace: `casp/services/bq_ops` → `scripts/bq_ops/entrypoints`
   - Verify workflow paths

3. **Test everything**
   - Run: `tox -e unit`
   - Run: `tox -e lint`
   - Fix any remaining import issues

4. **Update alembic.ini**
   - Verify migration script path references

### Important (Pre-merge)

5. **Update CLAUDE.md**
   - Document new structure
   - Update command examples
   - Update development setup

6. **Update README.rst**
   - Installation instructions with Poetry
   - New directory structure
   - Docker build instructions

7. **Remove old setup.py**
   - Poetry fully replaces setuptools
   - Keep only for legacy compatibility if needed

8. **Create migration guide**
   - For developers working on feature branches
   - How to rebase on new structure

### Nice to Have

9. **GitHub Actions / CI**
   - Update to use Poetry
   - Update Docker build workflows
   - Add `make requirements` step

10. **Documentation**
    - Update Sphinx conf if needed
    - Regenerate API docs

11. **Pre-commit hooks**
    - Add hook to regenerate requirements.txt on pyproject.toml change
    - Add import sorting with isort

## 📊 Statistics

- **Commits**: 4
  1. Add pyproject.toml + move deployments
  2. Restructure to cas_backend with apps/core
  3. Create namespace package + update imports
  4. Migrate to Poetry + update Docker

- **Files Modified**: 400+
- **Files Moved**: 218 (with `git mv` - history preserved)
- **Imports Updated**: 92 files automatically
- **Lines Changed**: ~10,000+

## 🏗️ New Structure

```
cellarium-cloud/repo/
├── cellarium/
│   └── cas_backend/
│       ├── apps/
│       │   ├── compute/          # Main REST API (CPU-intensive)
│       │   ├── admin/            # Flask admin UI
│       │   └── model_inference/  # ML model serving
│       └── core/
│           ├── auth/             # Authentication
│           ├── config.py         # Settings
│           ├── db/               # Database (PostgreSQL)
│           ├── data_managers/    # Data access (PG + BigQuery)
│           └── utils/            # Shared utilities
├── scripts/                      # Utility scripts
│   └── bq_ops/
│       ├── *.py                  # Main scripts
│       └── entrypoints/          # WDL workflow entry points
├── workflows/
│   ├── kubeflow/                 # ML pipelines
│   └── wdl/                      # Workflow definitions
├── deploy/
│   ├── docker/                   # Dockerfiles
│   ├── cloudrun/                 # GCP Cloud Run configs
│   └── scripts/                  # Deployment scripts
├── tests/
│   ├── unit/
│   └── integration/
├── pyproject.toml                # Poetry config (source of truth)
├── poetry.lock                   # Locked dependencies
├── requirements.txt              # Auto-generated (for Docker)
└── Makefile                      # Common commands
```

## 🔧 Commands Reference

### Development
```bash
# Install dependencies
make install  # or: poetry install

# Run tests
make test     # or: tox -e unit

# Lint code
make lint     # or: tox -e lint

# Format code
make format   # or: tox -e format
```

### Dependencies
```bash
# Add a new package
poetry add <package>

# Add dev dependency
poetry add --group dev <package>

# Update dependencies
poetry update

# Export for Docker
make requirements
```

### Docker
```bash
# Generate requirements.txt
make requirements

# Build images
docker build -f deploy/docker/Dockerfile.compute -t cas-compute .
docker build -f deploy/docker/Dockerfile.admin -t cas-admin .
docker build -f deploy/docker/Dockerfile.model_inference -t cas-ml .
```

### Database
```bash
# Run migrations
alembic -c cellarium/cas_backend/core/db/alembic.ini upgrade head

# Create migration
alembic -c cellarium/cas_backend/core/db/alembic.ini revision --autogenerate -m "description"
```

## 🚨 Breaking Changes

1. **All imports changed**
   - Old: `from casp.services.api import ...`
   - New: `from cellarium.cas_backend.apps.compute import ...`

2. **requirements.txt is auto-generated**
   - Don't edit manually
   - Run `make requirements` before Docker builds

3. **Directory structure completely changed**
   - `src/casp/` → `cellarium/cas_backend/`
   - `src/casp/services/api/` → `cellarium/cas_backend/apps/compute/`
   - `src/casp/services/db/` → `cellarium/cas_backend/core/db/`

4. **Docker builds changed**
   - New Dockerfiles in `deploy/docker/`
   - Must run `make requirements` first
   - Copy `cellarium/` not `src/`

5. **Settings location**
   - Still at `src/settings/.env` (temporary, consider moving)

## 📝 Notes

- Git history preserved for all moved files (used `git mv`)
- Namespace package allows future `cellarium.*` packages
- Poetry lock file ensures reproducible builds
- Old Dockerfiles kept temporarily for reference (can delete after verification)
