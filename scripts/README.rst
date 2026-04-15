Scripts
=======

Standalone utility scripts for one-off data preparation tasks. These are **not** part of the deployed
application and require the optional ``scripts`` dependency group::

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

Usage
~~~~~

.. code-block:: shell

    python scripts/create_cell_ontology_resource.py \
        --owl-url https://github.com/obophenotype/cell-ontology/releases/download/v2025-07-30/cl.owl \
        --output gs://your-bucket/path/to/cell_ontology_resource.json

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
