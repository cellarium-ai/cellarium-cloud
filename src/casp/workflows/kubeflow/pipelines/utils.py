import tempfile
import typing as t

from google.cloud import aiplatform
from kfp import compiler
from kfp.dsl.base_component import BaseComponent

from casp.services import utils
from casp.workflows.kubeflow import constants


def submit_pipeline(
    pipeline_func: t.Union[BaseComponent, t.Callable],
    pipeline_display_name: str,
    pipeline_kwargs: t.Dict[str, t.Any],
    pipeline_location: str = constants.DEFAULT_PIPELINE_LOCATION,
) -> None:
    """
    Create and run a pipeline on Vertex AI Pipelines. Use a temporary file to compile the pipeline config,
    then run the pipeline job and delete the temporary file.

    :param pipeline_func: Pipeline function, must be a function wrapped :func:`kfp.dsl.pipeline` decorator.
    :param pipeline_display_name: A name displayed in the Vertex AI Pipelines UI.
    :param pipeline_kwargs: Keyword arguments to pass to the pipeline function.
    :param pipeline_location: Datacenter location of Google Cloud Platform to run the pipeline job.
    """
    temp_file = tempfile.NamedTemporaryFile(suffix=".yaml")

    credentials, project_id = utils.get_google_service_credentials()
    aiplatform.init(project=project_id, location=pipeline_location, credentials=credentials)

    compiler.Compiler().compile(pipeline_func=pipeline_func, package_path=temp_file.name)

    job = aiplatform.PipelineJob(
        display_name=pipeline_display_name,
        template_path=temp_file.name,
        parameter_values=pipeline_kwargs,
    )

    job.submit(network=constants.NETWORK_NAME)

    temp_file.close()
