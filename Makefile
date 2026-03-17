.PHONY: help requirements docker-requirements clean install test lint format

help:
	@echo "Available commands:"
	@echo "  make requirements         - Export Poetry dependencies to requirements.txt for production"
	@echo "  make docker-requirements  - Export all dependencies (including dev) for Docker builds"
	@echo "  make install              - Install dependencies with Poetry"
	@echo "  make test                 - Run tests with tox"
	@echo "  make lint                 - Run linting"
	@echo "  make format               - Format code"
	@echo "  make clean                - Remove generated files"

# Export requirements.txt for production (Docker, etc.)
requirements:
	poetry export -f requirements.txt --output requirements.txt --without-hashes --without dev,test,docs

# Export requirements-dev.txt with all dependencies for development Docker builds
docker-requirements:
	poetry export -f requirements.txt --output requirements.txt --without-hashes --without dev,docs
	poetry export -f requirements.txt --output requirements-dev.txt --without-hashes --with dev,test

# Install dependencies
install:
	poetry install --with dev,test

# Run tests
test:
	poetry run tox -e unit

# Run linting
lint:
	poetry run tox -e lint

# Format code
format:
	poetry run tox -e format

# Clean generated files
clean:
	rm -rf .tox
	rm -rf htmlcov
	rm -rf .pytest_cache
	rm -rf .coverage coverage.xml
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	rm -f requirements.txt requirements-dev.txt
