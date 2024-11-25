import typer
from typing_extensions import Annotated

from casp.scripts.build_ontology_resources.benchmarking_resource import main as _build_benchmarking_resource
from casp.scripts.build_ontology_resources.consensus_engine_resource import main as _build_consensus_engine_resource

typer_app = typer.Typer()


@typer_app.command()
def build_consensus_engine_resource(
    cell_type_ontology_owl_file_url: Annotated[str, typer.Option()],
    output_file_path: Annotated[str, typer.Option()],
):
    """
    Build a consensus engine resource using the provided cell type ontology OWL file.

    This function processes the specified ontology file and outputs a resource file that can be used by
    the consensus engine for downstream tasks.

    :param cell_type_ontology_owl_file_url: URL to the cell type ontology OWL file
    :param output_file_path: Path to save the generated consensus engine resource

    Example usage
        To build a consensus engine resource, use the following command:

        .. code-block:: console

            python casp/scripts/build_consensus_engine_resource.py \\
                --cell-type-ontology-owl-file-url "https://example.com/ontology.owl" \\
                --output-file-path /path/to/output/resource.json
    """
    _build_consensus_engine_resource(
        cell_type_ontology_owl_file_url=cell_type_ontology_owl_file_url,
        output_file_path=output_file_path,
    )


@typer_app.command()
def build_benchmarking_resource(
    cell_type_ontology_owl_file_url: Annotated[str, typer.Option()],
    output_file_path: Annotated[str, typer.Option()],
    n_hops: Annotated[int, typer.Option()],
):
    """
    Create Cell Type Ontology Resources used in Ontology Aware Strategy and save it as a JSON file.

    :param cell_type_ontology_owl_file_url: URL to the cell type ontology OWL file
    :param output_file_path: Path to save the generated benchmarking resource
    :param n_hops: Number of hops for ontology-based neighborhood exploration

    Example usage
        To build a benchmarking resource, use the following command:

        .. code-block:: console

            python casp/scripts/build_benchmarking_resource.py \\
                --cell-type-ontology-owl-file-url "https://example.com/ontology.owl" \\
                --output-file-path /path/to/output/resource.json \\
                --n-hops 3
    """
    _build_benchmarking_resource(
        cell_type_ontology_owl_file_url=cell_type_ontology_owl_file_url,
        output_file_path=output_file_path,
        n_hops=n_hops,
    )


if __name__ == "__main__":
    typer_app()
