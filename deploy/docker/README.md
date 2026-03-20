# Docker Images

This directory contains Dockerfiles for building CAS Backend services.

## Images

- **Dockerfile.compute** - Main computational API service (formerly "api")
- **Dockerfile.admin** - Flask admin interface
- Cloud Run flavor settings live separately in `deploy/cloudrun/`

## Building Images

### Prerequisites

Generate `deploy/requirements.txt.lock` from `pyproject.toml` using Poetry:

```bash
# From repo root
make requirements
```

This creates `deploy/requirements.txt.lock` from Poetry dependencies (production only, no dev/test deps).

### Build Commands

```bash
# Compute service (main API)
docker build -f deploy/docker/Dockerfile.compute -t cas-backend-compute:latest .

# Admin service
docker build -f deploy/docker/Dockerfile.admin -t cas-backend-admin:latest .
```

### Image Structure

All images:
- Use Python 3.12
- Copy `cellarium/` package to `/app/cellarium/`
- Set `PYTHONPATH=/app` for module imports
- Run as non-root user `appuser` (UID 1000)
- Expect runtime configuration to be mounted at `/app/settings/.env`

## Running Containers

```bash
# Compute service (default port 8000)
docker run -p 8000:8000 \
  -v $(pwd)/settings/.env:/app/settings/.env \
  cas-backend-compute:latest

# Admin service (default port 5000)
docker run -p 5000:5000 \
  -v $(pwd)/settings/.env:/app/settings/.env \
  cas-backend-admin:latest
```

## CI/CD Integration

The build process should:

1. **Export requirements**: `make requirements`
2. **Build image**: `docker build -f deploy/docker/Dockerfile.<service> ...`
3. **Push to registry**: `docker push ...`

Example GitHub Actions workflow:

```yaml
- name: Export Poetry dependencies
  run: make requirements

- name: Build Docker image
  run: |
    docker build \
      -f deploy/docker/Dockerfile.compute \
      -t gcr.io/myproject/cas-backend-compute:${{ github.sha }} \
      .

- name: Push to GCR
  run: docker push gcr.io/myproject/cas-backend-compute:${{ github.sha }}
```

## Migration from Old Structure

Old Dockerfiles referenced:
- `COPY src .` → Now: `COPY cellarium/ ./cellarium/`
- `ENV PYTHONPATH=/app` → Now: `ENV PYTHONPATH=/app` (same, but imports use `cellarium.cas_backend`)
- `requirements.lock` manually maintained → Now: Generated via `make requirements` from `pyproject.toml` as `requirements.txt.lock`

## Troubleshooting

**Import errors**: Ensure `PYTHONPATH=/app` and imports use full path:
```python
# Correct
from cellarium.cas_backend.apps.compute import main

# Wrong
from cas_backend.apps.compute import main
```

**Missing dependencies**: Regenerate requirements.txt.lock:
```bash
make requirements
docker build --no-cache ...
```

**Settings not found**: Mount settings volume or copy into image:
```bash
docker run -v $(pwd)/settings/.env:/app/settings/.env ...
```
