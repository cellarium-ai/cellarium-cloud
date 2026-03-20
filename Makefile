.PHONY: help requirements docker-requirements clean install setup test lint format docs

help:
	@echo "Available commands:"
	@echo "  make requirements         - Export Poetry dependencies to deploy/requirements.txt.lock for production"
	@echo "  make docker-requirements  - Export production dependencies for Docker builds"
	@echo "  make setup                - Install Poetry"
	@echo "  make install              - Install dependencies with Poetry"
	@echo "  make test                 - Run unit tests"
	@echo "  make lint                 - Run linting"
	@echo "  make format               - Format code"
	@echo "  make docs                 - Build Sphinx documentation"
	@echo "  make clean                - Remove generated files"

# Export requirements.txt.lock for production (Docker, etc.)
requirements:
	poetry export -f requirements.txt --output deploy/requirements.txt.lock --without-hashes --without dev,test,docs

# Export requirements.txt.lock for Docker builds
docker-requirements:
	poetry export -f requirements.txt --output deploy/requirements.txt.lock --without-hashes --without dev,docs

# Install Poetry
setup:
	pip install --upgrade pip
	pip install poetry

# Install dependencies
install:
	poetry install --with dev,test

# Run tests
test:
	ENVIRONMENT=test poetry run pytest tests/unit

# Run linting
lint:
	poetry run ruff check --diff cellarium tests
	poetry run ruff format --diff cellarium tests

# Format code
format:
	poetry run ruff format cellarium tests
	poetry run ruff check --fix --unsafe-fixes cellarium tests

# Build documentation
docs:
	poetry run sphinx-build -W -b html docs/source docs/build/html

# Clean generated files
clean:
	rm -rf htmlcov
	rm -rf .pytest_cache
	rm -rf .coverage coverage.xml
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
