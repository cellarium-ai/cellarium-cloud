import typing as t

import typer
from typing_extensions import Annotated

from casp.scripts.benchmarking import calculate_metrics_for_cas_responses as _calculate_metrics
from casp.scripts.benchmarking import generate_cas_outputs as _generate_cas_outputs

typer_app = typer.Typer()


@typer_app.command()
def generate_cas_outputs(
    dataset_paths: Annotated[t.List[str], typer.Option()],
    model_name: Annotated[str, typer.Option()],
    cas_api_token: Annotated[str, typer.Option()],
    cas_api_url: Annotated[str, typer.Option()],
    cas_results_output_path: Annotated[str, typer.Option()],
):
    """
    Generate CAS responses for given datasets using the specified model and CAS API.

    Results will be saved as pickle files in `{cas_results_output_path}/{model_name}/` directory.

    :param dataset_paths: Paths to the dataset files. Can be a list of paths or a path to txt file with paths.
    :param model_name: Name of the model to use for generating CAS outputs
    :param cas_api_token: CAS API token for authentication
    :param cas_api_url: CAS API base URL
    :param cas_results_output_path: Path to save the CAS results output

    Example usage
        To generate CAS outputs, use the following command:

        .. code-block:: console

            python casp/scripts/generate_cas_outputs.py \\
                --dataset-paths "path/to/datasets.txt" \\
                --model-name my_model \\
                --cas-api-token my_api_token \\
                --cas-api-url "https://api.example.com" \\
                --cas-results-output-path /path/to/output/results.json
    """
    _generate_cas_outputs(
        dataset_paths=dataset_paths,
        model_name=model_name,
        cas_api_token=cas_api_token,
        cas_api_url=cas_api_url,
        cas_results_output_path=cas_results_output_path,
    )


@typer_app.command()
def calculate_metrics(
    dataset_paths: Annotated[t.Union[str, t.List[str]], typer.Option()],
    cas_result_paths: Annotated[t.Union[str, t.List[str]], typer.Option()],
    model_name: Annotated[str, typer.Option()],
    num_hops: Annotated[int, typer.Option()],
    co_resource_path: Annotated[str, typer.Option()],
    output_path: Annotated[str, typer.Option()],
    wandb_project: Annotated[str, typer.Option()],
    batch_size: Annotated[int, typer.Option()],
):
    """
    Calculate metrics for CAS responses and log results to Weights and Biases.

    :param dataset_paths: Paths to the dataset files. Can be a list of paths or a path to txt file with paths.
    :param cas_result_paths: Paths to the CAS result files. Can be a list of paths or a path to txt file with paths.
    :param model_name: Name of the model to evaluate
    :param num_hops: Number of hops for neighborhood exploration
    :param co_resource_path: Path to the consensus engine resource
    :param output_path: Path to save the calculated metrics
    :param wandb_project: Weights & Biases project name for logging results
    :param batch_size: Batch size for processing the data

    Example usage
        To calculate metrics, use the following command:

        .. code-block:: console

            python casp/scripts/calculate_metrics.py \\
                --dataset-paths "path/to/datasetpaths.txt" \\
                --cas-result-paths "path/to/results.txt"\\
                --model-name my_model \\
                --num-hops 3 \\
                --co-resource-path path/to/resource.json \\
                --output-path path/to/output/metrics.json \\
                --wandb-project my_project \\
                --batch-size 64
    """
    _calculate_metrics(
        dataset_paths=dataset_paths,
        cas_result_paths=cas_result_paths,
        model_name=model_name,
        num_hops=num_hops,
        co_resource_path=co_resource_path,
        output_path=output_path,
        wandb_project=wandb_project,
        batch_size=batch_size,
    )


if __name__ == "__main__":
    typer_app()
