Scripts
=======

Standalone utility scripts for one-off data preparation tasks. Each script has two
faces: a CLI entry point under ``scripts/`` and a library module importable from the
installed package under ``cellarium.cas_backend.scripts``.

Install the optional dependency group before running or importing::

    poetry install --with scripts

create_cell_ontology_resource.py
--------------------------------

Generates the Cell Type Ontology resource JSON consumed by the CAS API from a CL OWL file.

The output JSON contains five precomputed fields:

- ``ancestors_dictionary`` — maps each CL term to its ordered list of ancestors
- ``cell_ontology_term_id_to_cell_type`` — maps each CL term to its human-readable label
- ``children_dictionary`` — maps each CL term to its direct children
- ``shortest_path_lengths_from_cell_root`` — BFS distance from ``CL:0000000`` to each term
- ``longest_path_lengths_from_cell_root`` — longest DAG path from ``CL:0000000`` to each term

CLI usage
~~~~~~~~~

.. code-block:: shell

    poetry run create-cell-ontology-resource \
        --owl-url https://github.com/obophenotype/cell-ontology/releases/download/v2025-07-30/cl.owl \
        --output gs://your-bucket/path/to/cell_ontology_resource.json

Library import
~~~~~~~~~~~~~~

Requires the ``scripts`` extras (``owlready2``, ``networkx``):

.. code-block:: python

    from cellarium.cas_backend.scripts.create_cell_ontology_resource import (
        create_cell_ontology_resource,
    )

    create_cell_ontology_resource(
        owl_url="https://github.com/obophenotype/cell-ontology/releases/download/v2025-07-30/cl.owl",
        output="gs://your-bucket/path/to/cell_ontology_resource.json",
    )

Arguments
~~~~~~~~~

.. option:: --owl-url

   URL or local filesystem path to the CL OWL file. Accepts any URL supported by ``owlready2``.

.. option:: --output

   Destination file path for the resulting JSON. Accepts local paths or GCS URIs
   (``gs://bucket/path/resource.json``) via ``smart-open``.

Output
~~~~~~

The script writes a single JSON file. The resource must then be uploaded to the
``OntologicalColumn`` record in the database (via the Admin interface) and the GCS path stored
in ``OntologicalColumn.ontology_resource_name`` so the API can look it up by name.

create_vsindex.py
-----------------

Creates a TileDB IVF_FLAT vector search index from pre-computed .csv.gz embedding
batch files, in three phases: training ingest, streaming update, and consolidation.

CLI usage
~~~~~~~~~

.. code-block:: shell

    poetry run create-vsindex \
        --embeddings-prefix gs://bucket/path/transform \
        --index-path gs://bucket/path/ivflat_vsindex.soma \
        --total-batches 1234 \
        --embedding-dim 64

Library import
~~~~~~~~~~~~~~

Base-importable (only ``click`` is scripts-extras-only; the indexing deps are in main):

.. code-block:: python

    from cellarium.cas_backend.scripts.create_vsindex import create_vsindex

    create_vsindex(
        embeddings_prefix="gs://bucket/path/transform",
        index_path="gs://bucket/path/ivflat_vsindex.soma",
        total_batches=1234,
        embedding_dim=64,
    )

Arguments
~~~~~~~~~

.. option:: --embeddings-prefix

   Path prefix for ``batch_{i}.csv.gz`` files. Accepts local paths or GCS URIs
   (``gs://``). The script reads ``{prefix}/batch_0.csv.gz`` through
   ``{prefix}/batch_{N-1}.csv.gz``.

.. option:: --index-path

   Destination URI for the TileDB IVF_FLAT index group. Accepts local paths or
   ``gs://`` URIs.

.. option:: --total-batches

   Total number of batch files.

.. option:: --embedding-dim

   Dimensionality of each embedding vector (number of columns after the ID column).

.. option:: --training-sample-size

   Maximum number of vectors to load for the initial index creation. Default: 5 000 000.

.. option:: --distance-metric

   Distance metric: ``COSINE`` (default) or ``L2``.

.. option:: --normalize / --no-normalize

   Whether to L2-normalize embeddings before indexing. Default: enabled.

.. option:: --max-partitions

   Upper bound for IVF partition count (rounded to a power of two). Default: 65 536.

.. option:: --update-chunk-size

   Number of batches per streaming-update chunk. Default: 50.
