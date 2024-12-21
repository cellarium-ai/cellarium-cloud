"""
This module provides pytest fixtures for generating mock data used in unit tests.

The data includes cell metadata, mock AnnData objects, serialized AnnData files, and more. These fixtures are reusable
and help simulate different scenarios for testing functionalities that rely on biological or structured data.
"""

import io
import random
import tempfile
import typing as t

import anndata
import numpy as np
import pandas as pd
import pytest

from casp.services.api import schemas


@pytest.fixture
def cell_info_data() -> t.List[schemas.CellariumCellMetadata]:
    """
    Centralized fixture that represents the source of cell data. It is needed for adding cells into the database and to
    use the same cells for mocking MatchingEngine.

    :return: A list of instances of :class:`~casp.services.api.schemas.CellariumCellMetadata` with 'cell_type' and
        'cell_type_ontology_term_id'.
    """
    # Define a mapping between cell types and their ontology term IDs
    cell_type_to_ontology_id_map = {
        "monocyte": "CL:0000576",
        "erythrocyte": "CL:0000232",
        "lymphocyte": "CL:0000542",
        "CD4-positive, alpha-beta T cell": "CL:0000624",
    }

    return [
        schemas.CellariumCellMetadata(
            cas_cell_index=i, cell_type=cell_type, cell_type_ontology_term_id=cell_type_to_ontology_id_map[cell_type]
        )
        for i, cell_type in enumerate(random.choices(list(cell_type_to_ontology_id_map.keys()), k=40))
    ]


@pytest.fixture
def mock_valid_anndata() -> anndata.AnnData:
    """
    Fixture to provide a larger mocked AnnData object.

    :return: An `AnnData` object with random expression data and metadata.
    """
    n_cells = 20
    n_genes = 500

    # Generate obs dataframe (metadata for cells)
    obs_data = {
        "total_mrna_umis": np.random.uniform(50.0, 500.0, size=n_cells).astype(np.float32),  # Random float values
    }
    obs = pd.DataFrame(obs_data, index=[f"cell_{i}" for i in range(n_cells)])

    # Generate var dataframe (gene names)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    # Generate X matrix (gene expression data)
    X = np.random.rand(n_cells, n_genes).astype(np.float32)  # Random gene expression values

    return anndata.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def mock_file_with_anndata(mock_valid_anndata: anndata.AnnData) -> io.BytesIO:
    """
    Fixture to serialize the mocked AnnData object into a BytesIO object.

    :param mock_valid_anndata: The AnnData object to serialize.
    :return: A `BytesIO` object containing the serialized AnnData file.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temp_file:
        temp_file_path = temp_file.name
        # Write the AnnData object to a temporary file
        mock_valid_anndata.write_h5ad(temp_file_path)

    # Read the file back into a BytesIO object
    with open(temp_file_path, "rb") as f:
        file = io.BytesIO(f.read())

    return file
