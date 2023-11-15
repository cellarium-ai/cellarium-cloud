import typing as t

from kfp import dsl

from casp.workflows.vertex_ai.job_components.pca import create_pca_embed_job, create_pca_train_job
from casp.workflows.vertex_ai.job_components.registry import create_register_embedding_model_job


@dsl.pipeline(name="ipca_train", description="Incremental PCA Train")
def train_pipeline(gcs_config_path: str) -> None:
    """
    KFP pipeline to run PCA train pipeline.

    :param gcs_config_path: GCS path to the PCA train config file.
    """
    train_job_task_op = create_pca_train_job(gcs_config_path=gcs_config_path)
    _ = train_job_task_op()


@dsl.pipeline(name="ipca_embed", description="Incremental PCA Embed Data")
def embed_pipeline(gcs_config_path: str) -> None:
    """
    KFP pipeline to run PCA embed pipeline.

    :param gcs_config_path: GCS path to the PCA embed config file.
    """
    embed_job_task_op = create_pca_embed_job(gcs_config_path=gcs_config_path)
    _ = embed_job_task_op()


@dsl.pipeline(name="ipca_train_embed", description="Incremental PCA Train and Embed Data")
def train_embed_pipeline(train_gcs_config_path: str, embed_gcs_config_path: str) -> None:
    """
    KFP pipeline to run PCA train and embed pipelines in sequence.

    :param train_gcs_config_path: GCS path to the PCA train config file.
    :param embed_gcs_config_path: GCS path to the PCA embed config file.
    """
    train_job_task_op = create_pca_train_job(gcs_config_path=train_gcs_config_path)
    embed_job_task_op = create_pca_embed_job(gcs_config_path=embed_gcs_config_path)
    train_job_task = train_job_task_op()
    embed_job_task = embed_job_task_op()

    embed_job_task.after(train_job_task)


@dsl.pipeline(name="ipca_train_embed_bulk", description="Incremental PCA Train and Embed Data Bulk")
def train_embed_bulk_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run PCA train and embed pipelines in sequence for multiple models simultaneously.
    :class:`dsl.ParallelFor` doesn't work with a list of Python dictionaries, so we need to pass a list of JSON strings,
    which will be parsed into a list of dictionaries inside the pipeline.

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """
    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_embed_pipeline(
            train_gcs_config_path=item.train_gcs_config_path, embed_gcs_config_path=item.embed_gcs_config_path
        )
        # register_embedding_model_op = create_register_embedding_model_job(
        #     gcs_config_path=item.register_model_gcs_config_path
        # )
        # register_embedding_model_op()
