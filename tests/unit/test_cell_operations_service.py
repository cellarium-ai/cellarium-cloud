"""
Tests things in the cell_operations_service that can reasonably be tested with a unit test.
"""

import numpy

from casp.services.api.services.cell_operations_service import CellOperationsService

def test_split_embeddings_into_chunks_even_split() -> None:
    embeddings = numpy.random.rand(100, 100)
    chunk_size = 10

    # Turns out this is how you call a private method in a Python unit test.
    chunks = CellOperationsService._CellOperationsService__split_embeddings_into_chunks(embeddings, chunk_size)

    assert len(chunks) == 10
    assert sum(len(chunk) for chunk in chunks) == 100

def test_split_embeddings_into_chunks_uneven_split() -> None:
    embeddings = numpy.random.rand(100, 100)
    chunk_size = 9

    # Turns out this is how you call a private method in a Python unit test.
    chunks = CellOperationsService._CellOperationsService__split_embeddings_into_chunks(embeddings, chunk_size)

    assert len(chunks) == 12
    assert sum(len(chunk) for chunk in chunks) == 100