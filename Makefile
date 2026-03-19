.PHONY: help requirements docker-requirements clean install test lint format

help:
	@echo "Available commands:"
	@echo "  make requirements         - Export Poetry dependencies to deploy/requirements.txt.lock for production"
	@echo "  make docker-requirements  - Export all dependencies (including dev) for Docker builds"
	@echo "  make install              - Install dependencies with Poetry"
	@echo "  make test                 - Run unit tests"
	@echo "  make lint                 - Run linting"
	@echo "  make format               - Format code"
	@echo "  make clean                - Remove generated files"

# Export requirements.txt.lock for production (Docker, etc.)
requirements:
	poetry export -f requirements.txt --output deploy/requirements.txt.lock --without-hashes --without dev,test,docs

# Export requirements-dev.txt with all dependencies for development Docker builds
docker-requirements:
	poetry export -f requirements.txt --output deploy/requirements.txt.lock --without-hashes --without dev,docs
	poetry export -f requirements.txt --output deploy/requirements-dev.lock --without-hashes --with dev,test

# Install dependencies
install:
	poetry install --with dev,test

# Run tests
test:
	poetry run pytest tests/unit

# Run linting
lint:
	poetry run black --line-length 120 --check cellarium tests
	poetry run isort --check-only --diff --line-length 120 --profile black tests cellarium
	poetry run flake8 cellarium tests

# Format code
format:
	poetry run black --line-length 120 cellarium tests
	poetry run isort --line-length 120 --profile black tests cellarium

# Clean generated files
clean:
	rm -rf htmlcov
	rm -rf .pytest_cache
	rm -rf .coverage coverage.xml
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
