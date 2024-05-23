#!/usr/bin/env python

from setuptools import find_packages, setup

with open("requirements/base.txt") as fh:
    install_requires = fh.readlines()

with open("README.rst") as fh:
    long_description = fh.read()

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "0.1.0"
setup(
    name="cell-annotation-service-pilot",
    version=version,
    description="Cell Annotation Service Pilot",
    author="Variants Team",
    author_email="variants@broadinstitute.org",
    long_description=long_description,
    install_requires=install_requires,
    tests_require=["coverage", "pytest"],
    python_requires=">=3.10",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    include_package_data=True,
    extras_require={
        "docs": [
            "Pillow~=9.5",
            "Sphinx~=7.2",
            "sphinx_gallery~=0.14",
            "sphinx_rtd_theme~=1.3",
            "setuptools-git-versioning==1.13.5",
        ],
    },
)
