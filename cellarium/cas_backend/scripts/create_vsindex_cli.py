"""
Click CLI entry point for the TileDB IVF_FLAT vector search index creator.

Kept separate from create_vsindex.py so that the library module remains importable
without the ``scripts`` extras (click). Invoke via::

    poetry run create-vsindex --help
"""

import click

from cellarium.cas_backend.scripts.create_vsindex import create_vsindex


@click.command()
@click.option("--embeddings-prefix", required=True, help="Path prefix for batch_{i}.csv.gz files (local or gs://).")
@click.option("--index-path", required=True, help="Output TileDB index URI (local or gs://).")
@click.option(
    "--total-batches", required=True, type=int, help="Number of batch files (batch_0.csv.gz … batch_{N-1}.csv.gz)."
)
@click.option(
    "--embedding-dim", required=True, type=int, help="Vector dimensionality (number of embedding columns per row)."
)
@click.option(
    "--training-sample-size",
    default=5_000_000,
    show_default=True,
    type=int,
    help="Max vectors to load before creating the initial index.",
)
@click.option(
    "--distance-metric",
    default="COSINE",
    show_default=True,
    type=click.Choice(["COSINE", "L2"], case_sensitive=False),
    help="Distance metric for the IVF_FLAT index.",
)
@click.option(
    "--normalize/--no-normalize", default=True, show_default=True, help="L2-normalize embeddings before indexing."
)
@click.option(
    "--max-partitions", default=65536, show_default=True, type=int, help="Upper bound for IVF partition count."
)
@click.option(
    "--update-chunk-size", default=50, show_default=True, type=int, help="Number of batches per streaming-update chunk."
)
@click.option(
    "--allow-list-csv",
    default=None,
    type=str,
    help=(
        "Optional path (local or gs://) to a CSV with a single 'soma_joinid' column (header "
        "included). When provided, only cells whose ID appears in the allow-list are indexed; "
        "all others are dropped."
    ),
)
def main(
    embeddings_prefix: str,
    index_path: str,
    total_batches: int,
    embedding_dim: int,
    training_sample_size: int,
    distance_metric: str,
    normalize: bool,
    max_partitions: int,
    update_chunk_size: int,
    allow_list_csv: str | None,
) -> None:
    """Create a TileDB IVF_FLAT vector search index from .csv.gz embedding batches."""
    create_vsindex(
        embeddings_prefix=embeddings_prefix,
        index_path=index_path,
        total_batches=total_batches,
        embedding_dim=embedding_dim,
        training_sample_size=training_sample_size,
        distance_metric=distance_metric,
        normalize=normalize,
        max_partitions=max_partitions,
        update_chunk_size=update_chunk_size,
        allow_list_csv=allow_list_csv,
    )


if __name__ == "__main__":
    main()
