import typing as t

from kfp import dsl

from casp.workflows.vertex_ai.job_components.pca import create_pca_embed_job, create_pca_train_job
from casp.workflows.vertex_ai.job_components.registry import create_register_embedding_model_job
from casp.workflows.vertex_ai.job_components.space_vector_search_index import create_index_create_job


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

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """
    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_embed_pipeline(
            train_gcs_config_path=item.train_gcs_config_path, embed_gcs_config_path=item.embed_gcs_config_path
        )


@dsl.pipeline(name="ipca_full_cycle", description="Incremental PCA Full Cycle")
def model_full_cycle_pipeline(
    train_gcs_config_path: str,
    embed_gcs_config_path: str,
    register_model_gcs_config_path: str,
    index_create_config: str,
) -> None:
    """
    KFP pipeline to train a PCA model, embedd data, create and deploy index based on the model and register them in
    Cellarium Cloud database.

    :param train_gcs_config_path: GCS path to the PCA train config file.
    :param embed_gcs_config_path: GCS path to the PCA embed config file.
    :param register_model_gcs_config_path: GCS path to the PCA register model config file.
    :param index_create_config: GCS path to the PCA index create config file.
    """
    train_job_task_op = create_pca_train_job(gcs_config_path=train_gcs_config_path)
    embed_job_task_op = create_pca_embed_job(gcs_config_path=embed_gcs_config_path)
    register_embedding_model_op = create_register_embedding_model_job(gcs_config_path=register_model_gcs_config_path)
    index_create_op = create_index_create_job(gcs_config_path=index_create_config)

    train_job_task = train_job_task_op()
    embed_job_task = embed_job_task_op()
    register_embedding_model_task = register_embedding_model_op()
    index_create_task = index_create_op()

    embed_job_task.after(train_job_task)
    register_embedding_model_task.after(train_job_task)
    index_create_task.after(embed_job_task)
    index_create_task.after(register_embedding_model_task)


@dsl.pipeline(name="ipca_full_cycle_bulk", description="Incremental PCA Full Cycle Bulk")
def model_full_cycle_bulk_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to runa full cycle of Cellarium PCA model in sequence for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """
    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        model_full_cycle_pipeline(
            train_gcs_config_path=item.train_gcs_config_path,
            embed_gcs_config_path=item.embed_gcs_config_path,
            register_model_gcs_config_path=item.register_model_gcs_config_path,
            index_create_config=item.index_create_gcs_config_path,
        )
