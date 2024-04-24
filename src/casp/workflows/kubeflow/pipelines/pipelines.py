import typing as t

from kfp import dsl

from casp.workflows.kubeflow import constants, job_components
from casp.workflows.kubeflow.job_components.utils import create_job


@dsl.pipeline(name="pca_train", description="Incremental PCA Train")
def train_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run PCA train pipeline.

    :param pipeline_config_paths: List of JSON strings with train and embed config paths.
    """
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_job_task_op = create_job(
            dsl_component=job_components.pca.train,
            component_name=constants.PCA_TRAIN_COMPONENT_NAME,
            gcs_config_path=item.pca_train_gcs_config_path,
        )
        _ = train_job_task_op()


@dsl.pipeline(name="pca_embed_parallel", description="Incremental PCA Embed Data")
def embed_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run PCA embed pipeline.

    :param pipeline_config_paths: List of JSON strings with train and embed config paths.
    """
    with dsl.ParallelFor(pipeline_config_paths) as item:
        embed_job_task_op = create_job(
            dsl_component=job_components.pca.embed,
            component_name=constants.PCA_EMBED_COMPONENT_NAME,
            gcs_config_path=item.pca_embed_gcs_config_path,
        )
        _ = embed_job_task_op()


@dsl.pipeline(name="pca_train_embed_parallel", description="Incremental PCA Train and Embed Data Parallel")
def train_embed_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run PCA train and embed pipelines in sequence for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """

    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_job_task_op = create_job(
            dsl_component=job_components.pca.train,
            component_name=constants.PCA_TRAIN_COMPONENT_NAME,
            gcs_config_path=item.pca_train_gcs_config_path,
        )
        embed_job_task_op = create_job(
            dsl_component=job_components.pca.embed,
            component_name=constants.PCA_EMBED_COMPONENT_NAME,
            gcs_config_path=item.pca_embed_gcs_config_path,
        )

        train_job_task = train_job_task_op()
        embed_job_task = embed_job_task_op()

        embed_job_task.after(train_job_task)


@dsl.pipeline(name="pca_deploy_index_parallel", description="PCA Deploy Space Vector Index Parallel")
def pca_deploy_index_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run PCA deploy index pipeline for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """

    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        index_create_op = create_job(
            dsl_component=job_components.pca_index_create.create_deploy_register_index,
            component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
            gcs_config_path=item.pca_index_create_gcs_config_path,
        )
        _ = index_create_op()


@dsl.pipeline(name="pca_full_cycle_parallel", description="Incremental PCA Full Cycle Parallel")
def pca_full_cycle_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to runa full cycle of Cellarium PCA model in sequence for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
    """

    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_job_task_op = create_job(
            dsl_component=job_components.pca.train,
            component_name=constants.PCA_TRAIN_COMPONENT_NAME,
            gcs_config_path=item.pca_train_gcs_config_path,
        )
        embed_job_task_op = create_job(
            dsl_component=job_components.pca.embed,
            component_name=constants.PCA_EMBED_COMPONENT_NAME,
            gcs_config_path=item.pca_embed_gcs_config_path,
        )
        register_embedding_model_op = create_job(
            dsl_component=job_components.registry.register_embedding_model,
            component_name=constants.PCA_REGISTRY_COMPONENT_NAME,
            gcs_config_path=item.pca_register_model_gcs_config_path,
        )
        index_create_op = create_job(
            dsl_component=job_components.pca_index_create.create_deploy_register_index,
            component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
            gcs_config_path=item.pca_index_create_gcs_config_path,
        )

        train_job_task = train_job_task_op()
        embed_job_task = embed_job_task_op()
        register_embedding_model_task = register_embedding_model_op()
        index_create_task = index_create_op()

        embed_job_task.after(train_job_task)
        register_embedding_model_task.after(train_job_task)
        index_create_task.after(embed_job_task)
        index_create_task.after(register_embedding_model_task)


@dsl.pipeline(name="pca_resize_full_cycle_parallel", description="Incremental PCA Resize and Save Full Cycle Parallel")
def pca_resize_full_cycle_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline that takes a model resizes it, saves and the uses for creating new index, and then registers it.

    :param pipeline_config_paths: List of JSON strings with train and embed config paths.
    """
    with dsl.ParallelFor(pipeline_config_paths) as item:
        resize_and_save_op = create_job(
            dsl_component=job_components.pca.resize_and_save,
            component_name=constants.PCA_RESIZE_AND_SAVE_COMPONENT_NAME,
            gcs_config_path=item.pca_resize_and_save_gcs_config_path,
        )
        embed_job_task_op = create_job(
            dsl_component=job_components.pca.embed,
            component_name=constants.PCA_EMBED_COMPONENT_NAME,
            gcs_config_path=item.pca_embed_gcs_config_path,
        )
        register_embedding_model_op = create_job(
            dsl_component=job_components.registry.register_embedding_model,
            component_name=constants.PCA_REGISTRY_COMPONENT_NAME,
            gcs_config_path=item.pca_register_model_gcs_config_path,
        )
        index_create_op = create_job(
            dsl_component=job_components.pca_index_create.create_deploy_register_index,
            component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
            gcs_config_path=item.pca_index_create_gcs_config_path,
        )
        resize_and_save_task = resize_and_save_op()
        embed_job_task = embed_job_task_op()
        register_embedding_model_task = register_embedding_model_op()
        index_create_task = index_create_op()

        embed_job_task.after(resize_and_save_task)
        register_embedding_model_task.after(resize_and_save_task)
        index_create_task.after(embed_job_task)
        index_create_task.after(register_embedding_model_task)


@dsl.pipeline(name="summary_stats_train_parallel", description="Summary Stats Train in Parallel")
def summary_stats_train_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run Summary Stats train pipeline for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train config paths.
    """

    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_mean_var_std = create_job(
            dsl_component=job_components.summary_stats.mean_var_std_train,
            component_name=constants.MEAN_VAR_STD_COMPONENT_NAME,
            gcs_config_path=item.mean_var_std_train_gcs_config_path,
        )
        train_tdigest = create_job(
            dsl_component=job_components.summary_stats.tdigest_train,
            component_name=constants.TDIGEST_COMPONENT_NAME,
            gcs_config_path=item.tdigest_train_gcs_config_path,
        )
        create_filter_feature_for_tdigest_op = create_job(
            dsl_component=job_components.summary_stats.tdigest_filter_features_component,
            component_name=constants.TDIGEST_FILTER_FEATURES_COMPONENT_NAME,
            gcs_config_path=item.tdigest_filter_features_gcs_config_path,
        )

        _ = train_mean_var_std()
        train_tdigest_task = train_tdigest()
        create_filter_feature_for_tdigest_task = create_filter_feature_for_tdigest_op()
        create_filter_feature_for_tdigest_task.after(train_tdigest_task)


@dsl.pipeline(name="logistic_regression_train_parallel", description="Logistic Regression Train Parallel")
def logistic_regression_train_pipeline(pipeline_config_paths: t.List[str]) -> None:
    """
    KFP pipeline to run Logistic Regression train pipeline for multiple models simultaneously.

    :param pipeline_config_paths:  List of JSON strings with train config paths.
    """

    with dsl.ParallelFor(pipeline_config_paths) as item:
        train_logistic_regression = create_job(
            dsl_component=job_components.logistic_regression.train_component,
            component_name=constants.LOGISTIC_REGRESSION_TRAIN_COMPONENT_NAME,
            gcs_config_path=item.logistic_regression_train_gcs_config_path,
        )

        _ = train_logistic_regression()


@dsl.pipeline(
    name="pca_train_base_and_resize_full_cycle_parallel", description="PCA Train Base and Resize Full Cycle Parallel"
)
def pca_train_base_and_resize_full_cycle_pipeline(
    train_base_config_paths: t.List[str], resize_full_cycle_config_paths: t.List[str]
) -> None:
    train_pipeline_task = train_pipeline(pipeline_config_paths=train_base_config_paths)
    resize_full_cycle_task = pca_resize_full_cycle_pipeline(pipeline_config_paths=resize_full_cycle_config_paths)
    resize_full_cycle_task.after(train_pipeline_task)


@dsl.pipeline(
    name="summary_stats_lr_train_parallel", description="Summary Stats and Logistic Regression Train Parallel"
)
def summary_stats_lr_train_pipeline(
    summary_stats_config_paths: t.List[str], logistic_regression_config_paths: t.List[str]
) -> None:
    summary_stats_train_task = summary_stats_train_pipeline(pipeline_config_paths=summary_stats_config_paths)
    logistic_regression_train_task = logistic_regression_train_pipeline(
        pipeline_config_paths=logistic_regression_config_paths
    )
    logistic_regression_train_task.after(summary_stats_train_task)


@dsl.pipeline(name="benchmark_cas_parallel", description="Benchmarking CAS Parallel")
def benchmark_cas_pipeline(pipeline_config_paths: t.List[str]) -> None:
    with dsl.ParallelFor(pipeline_config_paths) as item:
        benchmark_cas_op = create_job(
            dsl_component=job_components.benchmarking.benchmark_cas,
            component_name=constants.BENCHMARKING_COMPONENT_NAME,
            gcs_config_path=item.benchmarking_gcs_config_path,
        )
        _ = benchmark_cas_op()