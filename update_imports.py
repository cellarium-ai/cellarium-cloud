#!/usr/bin/env python3
"""
Script to update all imports from 'casp' to 'cellarium.cas_backend' with new structure.

This handles the major refactoring:
- casp.services.api → cellarium.cas_backend.apps.compute
- casp.services.admin → cellarium.cas_backend.apps.admin
- casp.services.model_inference → cellarium.cas_backend.apps.model_inference
- casp.services.db → cellarium.cas_backend.core.db
- casp.services._auth → cellarium.cas_backend.core.auth
- casp.services._settings → cellarium.cas_backend.core.config
- casp.services → cellarium.cas_backend.core
- casp.data_manager → cellarium.cas_backend.core.data_managers
- casp.scripts → scripts (now at repo root, not a package)
"""

import os
import re
from pathlib import Path


def update_imports_in_file(file_path):
    """Update all casp imports in a single file."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    original_content = content

    # Order matters! More specific patterns first
    replacements = [
        # Apps
        (r'from casp\.services\.api\.', 'from cellarium.cas_backend.apps.compute.'),
        (r'import casp\.services\.api\.', 'import cellarium.cas_backend.apps.compute.'),
        (r'from casp\.services\.admin\.', 'from cellarium.cas_backend.apps.admin.'),
        (r'import casp\.services\.admin\.', 'import cellarium.cas_backend.apps.admin.'),
        (r'from casp\.services\.model_inference\.', 'from cellarium.cas_backend.apps.model_inference.'),
        (r'import casp\.services\.model_inference\.', 'import cellarium.cas_backend.apps.model_inference.'),

        # Core - specific modules first
        (r'from casp\.services\.db\.', 'from cellarium.cas_backend.core.db.'),
        (r'import casp\.services\.db\.', 'import cellarium.cas_backend.core.db.'),
        (r'from casp\.services\.db import', 'from cellarium.cas_backend.core.db import'),
        (r'import casp\.services\.db$', 'import cellarium.cas_backend.core.db'),

        (r'from casp\.services\._auth\.', 'from cellarium.cas_backend.core.auth.'),
        (r'import casp\.services\._auth\.', 'import cellarium.cas_backend.core.auth.'),
        (r'from casp\.services\._auth import', 'from cellarium.cas_backend.core.auth import'),

        (r'from casp\.services import settings', 'from cellarium.cas_backend.core import config as settings'),
        (r'from casp\.services\._settings', 'from cellarium.cas_backend.core.config'),
        (r'import casp\.services\._settings', 'import cellarium.cas_backend.core.config'),

        (r'from casp\.services\.utils\.', 'from cellarium.cas_backend.core.utils.'),
        (r'import casp\.services\.utils\.', 'import cellarium.cas_backend.core.utils.'),

        # Generic services (after specific ones)
        (r'from casp\.services\.', 'from cellarium.cas_backend.core.'),
        (r'import casp\.services\.', 'import cellarium.cas_backend.core.'),
        (r'from casp\.services import', 'from cellarium.cas_backend.core import'),
        (r'import casp\.services', 'import cellarium.cas_backend.core'),

        # Data managers
        (r'from casp\.services\.api\.data_manager\.', 'from cellarium.cas_backend.core.data_managers.'),
        (r'import casp\.services\.api\.data_manager\.', 'import cellarium.cas_backend.core.data_managers.'),
        (r'from casp\.data_manager\.', 'from cellarium.cas_backend.core.data_managers.'),
        (r'import casp\.data_manager\.', 'import cellarium.cas_backend.core.data_managers.'),
        (r'from casp\.data_manager import', 'from cellarium.cas_backend.core.data_managers import'),

        # Scripts - now top-level, not a package import
        (r'from casp\.scripts\.', 'from scripts.'),
        (r'import casp\.scripts\.', 'import scripts.'),

        # Generic casp catch-all (should be last)
        (r'from casp\.', 'from cellarium.cas_backend.'),
        (r'import casp\.', 'import cellarium.cas_backend.'),
        (r'\bcasp\.', 'cellarium.cas_backend.'),  # For string references
    ]

    for pattern, replacement in replacements:
        content = re.sub(pattern, replacement, content)

    # Special case: settings module name changed
    content = content.replace('from cellarium.cas_backend.core import config as settings',
                              'from cellarium.cas_backend.core.config import settings')
    content = content.replace('settings.API_SERVICE_PORT', 'settings.API_SERVICE_PORT')

    if content != original_content:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False


def main():
    """Update all Python files in the project."""
    repo_root = Path(__file__).parent
    python_files = []

    # Find all Python files in relevant directories
    for directory in ['cellarium', 'scripts', 'tests', 'workflows']:
        dir_path = repo_root / directory
        if dir_path.exists():
            python_files.extend(dir_path.rglob('*.py'))

    # Also check setup.py, tox.ini, etc.
    for file in ['setup.py']:
        file_path = repo_root / file
        if file_path.exists():
            python_files.append(file_path)

    updated_files = []
    for file_path in python_files:
        if '__pycache__' in str(file_path) or '.tox' in str(file_path):
            continue

        if update_imports_in_file(file_path):
            updated_files.append(file_path)
            print(f"✓ Updated: {file_path.relative_to(repo_root)}")

    print(f"\n{'='*60}")
    print(f"Total files updated: {len(updated_files)}")
    print(f"{'='*60}")

    if updated_files:
        print("\nNext steps:")
        print("1. Review changes: git diff")
        print("2. Run tests: tox -e unit")
        print("3. Fix any remaining import errors manually")
        print("4. Commit: git commit -am 'refactor: update imports to cellarium.cas_backend'")


if __name__ == '__main__':
    main()
