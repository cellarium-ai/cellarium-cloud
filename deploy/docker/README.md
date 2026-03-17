# Docker Images

This directory contains Dockerfiles for building CAS Backend services.

## Images

- **Dockerfile.compute** - Main computational API service (formerly "api")
- **Dockerfile.admin** - Flask admin interface
- **Dockerfile.model_inference** - ML model serving
- **Dockerfile.pytorch** - Legacy PyTorch base (deprecated, use Dockerfile.model_inference)
- **Dockerfile.pytorch_cuda** - Legacy PyTorch with CUDA (deprecated)

## Building Images

### Prerequisites

Generate `requirements.txt` from `pyproject.toml` using Poetry:

```bash
# From repo root
make requirements
```

This creates `requirements.txt` from Poetry dependencies (production only, no dev/test deps).

### Build Commands

```bash
# Compute service (main API)
docker build -f deploy/docker/Dockerfile.compute -t cas-backend-compute:latest .

# Admin service
docker build -f deploy/docker/Dockerfile.admin -t cas-backend-admin:latest .

# Model inference service
docker build -f deploy/docker/Dockerfile.model_inference -t cas-backend-model-inference:latest .
```

### Image Structure

All images:
- Use Python 3.11
- Copy `cellarium/` package to `/app/cellarium/`
- Copy `scripts/` to `/app/scripts/` (compute service only)
- Copy `src/settings/` for environment configuration
- Set `PYTHONPATH=/app` for module imports
- Run as non-root user `appuser` (UID 1000)

## Running Containers

```bash
# Compute service (default port 8000)
docker run -p 8000:8000 \
  -v $(pwd)/src/settings/.env:/app/src/settings/.env \
  cas-backend-compute:latest

# Admin service (default port 5000)
docker run -p 5000:5000 \
  -v $(pwd)/src/settings/.env:/app/src/settings/.env \
  cas-backend-admin:latest

# Model inference service
docker run -p 8001:8001 \
  -v $(pwd)/src/settings/.env:/app/src/settings/.env \
  --gpus all \  # If using CUDA
  cas-backend-model-inference:latest
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
- `requirements.txt` manually maintained → Now: Generated via `make requirements` from `pyproject.toml`

## Troubleshooting

**Import errors**: Ensure `PYTHONPATH=/app` and imports use full path:
```python
# Correct
from cellarium.cas_backend.apps.compute import main

# Wrong
from cas_backend.apps.compute import main
```

**Missing dependencies**: Regenerate requirements.txt:
```bash
make requirements
docker build --no-cache ...
```

**Settings not found**: Mount settings volume or copy into image:
```bash
docker run -v $(pwd)/src/settings/.env:/app/src/settings/.env ...
```
