import typing as t

from kfp import dsl

from casp.workflows.kubeflow import constants, job_components
from casp.workflows.kubeflow.job_components.utils import create_job
from casp.workflows.kubeflow.pipelines import utils


@dsl.pipeline(name="pca_train", description="Incremental PCA Train")
def pca_train_pipeline(pca_train_pipeline_config_paths: t.List[str]):
    """
    KFP pipeline to run PCA train pipeline.

    :param pca_train_pipeline_config_paths: List of JSON strings with train and embed config paths.

    :rtype: dsl.Pipeline
    """
    with dsl.ParallelFor(pca_train_pipeline_config_paths, parallelism=54) as item:
        train_job_task_op = create_job(
            component_func=job_components.pca.train,
            component_name=constants.PCA_TRAIN_COMPONENT_NAME,
            gcs_config_path=item.pca_train_gcs_config_path,
        )
        _ = train_job_task_op()


# @dsl.pipeline(name="pca_embed_parallel", description="Incremental PCA Embed Data")
# def embed_pipeline(pca_embed_pipeline_config_paths: t.List[str]) -> None:
#     """
#     KFP pipeline to run PCA embed pipeline.
#
#     :param pca_embed_pipeline_config_paths: List of JSON strings with train and embed config paths.
#     """
#     with dsl.ParallelFor(pca_embed_pipeline_config_paths) as item:
#         embed_job_task_op = create_job(
#             component_func=job_components.pca.embed,
#             component_name=constants.PCA_EMBED_COMPONENT_NAME,
#             gcs_config_path=item.pca_embed_gcs_config_path,
#         )
#         _ = embed_job_task_op()
#
#
# @dsl.pipeline(name="pca_train_embed_parallel", description="Incremental PCA Train and Embed Data Parallel")
# def train_embed_pipeline(pipeline_config_paths: t.List[str]) -> None:
#     """
#     KFP pipeline to run PCA train and embed pipelines in sequence for multiple models simultaneously.
#
#     :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
#     """
#
#     # :class:`dsl.ParallelFor` requires a list of JSON strings as input
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         train_job_task_op = create_job(
#             component_func=job_components.pca.train,
#             component_name=constants.PCA_TRAIN_COMPONENT_NAME,
#             gcs_config_path=item.pca_train_gcs_config_path,
#         )
#         embed_job_task_op = create_job(
#             component_func=job_components.pca.embed,
#             component_name=constants.PCA_EMBED_COMPONENT_NAME,
#             gcs_config_path=item.pca_embed_gcs_config_path,
#         )
#
#         train_job_task = train_job_task_op()
#         embed_job_task = embed_job_task_op()
#
#         embed_job_task.after(train_job_task)
#
#
# @dsl.pipeline(name="pca_deploy_index_parallel", description="PCA Deploy Space Vector Index Parallel")
# def pca_deploy_index_pipeline(pipeline_config_paths: t.List[str]) -> None:
#     """
#     KFP pipeline to run PCA deploy index pipeline for multiple models simultaneously.
#
#     :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
#     """
#
#     # :class:`dsl.ParallelFor` requires a list of JSON strings as input
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         index_create_op = create_job(
#             component_func=job_components.pca_index_create.create_register_index,
#             component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
#             gcs_config_path=item.pca_index_create_gcs_config_path,
#         )
#         _ = index_create_op()
#
####
# @dsl.pipeline(name="pca_full_cycle_parallel", description="Incremental PCA Full Cycle Parallel")
# def pca_full_cycle_pipeline(pipeline_config_paths: t.List[str]) -> None:
#     """
#     KFP pipeline to runa full cycle of Cellarium PCA model in sequence for multiple models simultaneously.
#
#     :param pipeline_config_paths:  List of JSON strings with train and embed config paths.
#     """
#
#     # :class:`dsl.ParallelFor` requires a list of JSON strings as input
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         train_job_task_op = create_job(
#             component_func=job_components.pca.train,
#             component_name=constants.PCA_TRAIN_COMPONENT_NAME,
#             gcs_config_path=item.pca_train_gcs_config_path,
#         )
#         embed_job_task_op = create_job(
#             component_func=job_components.pca.embed,
#             component_name=constants.PCA_EMBED_COMPONENT_NAME,
#             gcs_config_path=item.pca_embed_gcs_config_path,
#         )
#         register_embedding_model_op = create_job(
#             component_func=job_components.registry.register_embedding_model,
#             component_name=constants.PCA_REGISTRY_COMPONENT_NAME,
#             gcs_config_path=item.pca_register_model_gcs_config_path,
#         )
#         index_create_op = create_job(
#             component_func=job_components.pca_index_create.create_register_index,
#             component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
#             gcs_config_path=item.pca_index_create_gcs_config_path,
#         )
#
#         train_job_task = train_job_task_op()
#         embed_job_task = embed_job_task_op()
#         register_embedding_model_task = register_embedding_model_op()
#         index_create_task = index_create_op()
#
#         embed_job_task.after(train_job_task)
#         register_embedding_model_task.after(train_job_task)
#         index_create_task.after(embed_job_task)
#         index_create_task.after(register_embedding_model_task)


# @dsl.pipeline(name="pca_resize_and_embed_parallel", description="Incremental PCA Resize and Embed Parallel")
# def pca_resize_and_embed_pipeline(pca_resize_and_embed_pipeline_config_paths: t.List[str]):
#     """
#     KFP pipeline that takes a model resizes it, saves and the uses for creating new index, and then registers it.
#
#     :param pca_resize_and_embed_pipeline_config_paths: List of JSON strings with train and embed config paths.
#
#     :rtype: dsl.Pipeline
#     """
#     with dsl.ParallelFor(pca_resize_and_embed_pipeline_config_paths, parallelism=72) as item:
#         resize_and_save_op = create_job(
#             component_func=job_components.pca.resize_and_save,
#             component_name=constants.PCA_RESIZE_AND_SAVE_COMPONENT_NAME,
#             gcs_config_path=item.pca_resize_and_save_gcs_config_path,
#         )
#         # embed_job_task_op = create_job(
#         #     component_func=job_components.pca.embed,
#         #     component_name=constants.PCA_EMBED_COMPONENT_NAME,
#         #     gcs_config_path=item.pca_embed_gcs_config_path,
#         # )
#         resize_and_save_task = resize_and_save_op()
#         # embed_job_task = embed_job_task_op()
#         # embed_job_task.after(resize_and_save_task)
@dsl.pipeline(name="pca_resize_in_parallel", description="Incremental PCA Resize ")
def pca_resize_pipeline(pca_resize_pipeline_config_paths: t.List[str]):
    """
    KFP pipeline that takes a model resizes it, saves and the uses for creating new index, and then registers it.

    :param pca_resize_pipeline_config_paths: List of JSON strings with train and embed config paths.

    :rtype: dsl.Pipeline
    """
    with dsl.ParallelFor(pca_resize_pipeline_config_paths) as item:
        resize_and_save_op = create_job(
            component_func=job_components.pca.resize_and_save,
            component_name=constants.PCA_RESIZE_AND_SAVE_COMPONENT_NAME,
            gcs_config_path=item.pca_resize_and_save_gcs_config_path,
        )
        _ = resize_and_save_op()

# def pca_resize_full_cycle_pipeline(pca_resize_full_cycle_pipeline_config_paths: t.List[str]):
#     """
#     KFP pipeline that takes a model resizes it, saves and the uses for creating new index, and then registers it.
#
#     :param pca_resize_full_cycle_pipeline_config_paths: List of JSON strings with train and embed config paths.
#
#     :rtype: dsl.Pipeline
#     """
#     with dsl.ParallelFor(pca_resize_full_cycle_pipeline_config_paths, parallelism=3) as item:
#         resize_and_save_op = create_job(
#             component_func=job_components.pca.resize_and_save,
#             component_name=constants.PCA_RESIZE_AND_SAVE_COMPONENT_NAME,
#             gcs_config_path=item.pca_resize_and_save_gcs_config_path,
#         )
        # embed_job_task_op = create_job(
        #     component_func=job_components.pca.embed,
        #     component_name=constants.PCA_EMBED_COMPONENT_NAME,
        #     gcs_config_path=item.pca_embed_gcs_config_path,
        # )
        # register_embedding_model_op = create_job(
        #     component_func=job_components.registry.register_embedding_model,
        #     component_name=constants.PCA_REGISTRY_COMPONENT_NAME,
        #     gcs_config_path=item.pca_registry_gcs_config_path,
        # )
        # index_create_op = create_job(
        #     component_func=job_components.pca_index_create.create_register_index,
        #     component_name=constants.PCA_INDEX_CREATE_COMPONENT_NAME,
        #     gcs_config_path=item.pca_index_create_gcs_config_path,
        # )
        # resize_and_save_task = resize_and_save_op()
        # embed_job_task = embed_job_task_op()
        # register_embedding_model_task = register_embedding_model_op()
        # index_create_task = index_create_op()

        # embed_job_task.after(resize_and_save_task)
        # register_embedding_model_task.after(resize_and_save_task)
        # index_create_task.after(embed_job_task)
        # index_create_task.after(register_embedding_model_task)


@dsl.pipeline(name="summary_stats_train_parallel", description="Summary Stats Train in Parallel")
def summary_stats_train_pipeline(summary_stats_train_pipeline_config_paths: t.List[str]):
    """
    KFP pipeline to run Summary Stats train pipeline for multiple models simultaneously.

    :param summary_stats_train_pipeline_config_paths:  List of JSON strings with train config paths.

    :rtype: dls.Pipeline
    """

    # :class:`dsl.ParallelFor` requires a list of JSON strings as input
    with dsl.ParallelFor(summary_stats_train_pipeline_config_paths) as item:
        train_mean_var_std = create_job(
            component_func=job_components.summary_stats.mean_var_std_train,
            component_name=constants.MEAN_VAR_STD_COMPONENT_NAME,
            gcs_config_path=item.mean_var_std_train_gcs_config_path,
        )
        train_tdigest = create_job(
            component_func=job_components.summary_stats.tdigest_train,
            component_name=constants.TDIGEST_COMPONENT_NAME,
            gcs_config_path=item.tdigest_train_gcs_config_path,
        )
        create_filter_feature_for_tdigest_op = create_job(
            component_func=job_components.summary_stats.tdigest_filter_features_component,
            component_name=constants.TDIGEST_FILTER_FEATURES_COMPONENT_NAME,
            gcs_config_path=item.tdigest_filter_features_gcs_config_path,
        )

        _ = train_mean_var_std()
        train_tdigest_task = train_tdigest()
        create_filter_feature_for_tdigest_task = create_filter_feature_for_tdigest_op()
        create_filter_feature_for_tdigest_task.after(train_tdigest_task)


#
#
# @dsl.pipeline(name="logistic_regression_train_parallel", description="Logistic Regression Train Parallel")
# def logistic_regression_train_pipeline(logistic_regression_config_paths: t.List[str]):
#     """
#     KFP pipeline to run Logistic Regression train pipeline for multiple models simultaneously.
#
#     :param logistic_regression_config_paths:  List of JSON strings with train config paths.
#     :rtype: dls.Pipeline
#     """
#
#     with dsl.ParallelFor(logistic_regression_config_paths) as item:
#         train_logistic_regression = create_job(
#             component_func=job_components.logistic_regression.train_component,
#             component_name=constants.LOGISTIC_REGRESSION_TRAIN_COMPONENT_NAME,
#             gcs_config_path=item.logistic_regression_train_gcs_config_path,
#         )
#
#         _ = train_logistic_regression()
#
@dsl.pipeline(
    name="pca_train_base_and_resize_parallel", description="PCA Train Base and Resize Full Cycle Parallel"
)
def pca_train_base_and_resize_pipeline(
        pca_train_pipeline_config_paths: t.List[str], pca_resize_pipeline_config_paths: t.List[str]
):
    train_pipeline_task = pca_train_pipeline(pca_train_pipeline_config_paths=pca_train_pipeline_config_paths)
    resize_full_task = pca_resize_pipeline(
        pca_resize_pipeline_config_paths=pca_resize_pipeline_config_paths
    )
    resize_full_task.after(train_pipeline_task)

# @dsl.pipeline(
#     name="pca_train_base_and_resize_full_cycle_parallel", description="PCA Train Base and Resize Full Cycle Parallel"
# )
# def pca_train_base_and_resize_full_cycle_pipeline(
#         pca_train_pipeline_config_paths: t.List[str], pca_resize_full_cycle_pipeline_config_paths: t.List[str]
# ):
#     train_pipeline_task = pca_train_pipeline(pca_train_pipeline_config_paths=pca_train_pipeline_config_paths)
#     resize_full_cycle_task = pca_resize_full_cycle_pipeline(
#         pca_resize_full_cycle_pipeline_config_paths=pca_resize_full_cycle_pipeline_config_paths
#     )
#     resize_full_cycle_task.after(train_pipeline_task)


#
#
# @dsl.pipeline(
#     name="summary_stats_lr_train_parallel", description="Summary Stats and Logistic Regression Train Parallel"
# )
# def summary_stats_lr_train_pipeline(
#     summary_stats_config_paths: t.List[str], logistic_regression_config_paths: t.List[str]
# ) -> None:
#     summary_stats_train_task = summary_stats_train_pipeline(pipeline_config_paths=summary_stats_config_paths)
#     logistic_regression_train_task = logistic_regression_train_pipeline(
#         logistic_regression_config_paths=logistic_regression_config_paths
#     )
#     logistic_regression_train_task.after(summary_stats_train_task)
#
#
# @dsl.pipeline(name="generate_cas_outputs", description="Generate outputs from CAS by running on the set of datasets")
def generate_cas_outputs_pipeline(generate_cas_outputs_pipeline_config_paths: t.List[str]):
    with dsl.ParallelFor(generate_cas_outputs_pipeline_config_paths, parallelism=18) as item:
        generate_cas_outputs_op = create_job(
            component_func=job_components.benchmarking.generate_cas_outputs,
            component_name=constants.GENERATE_CAS_OUTPUTS_COMPONENT_NAME,
            gcs_config_path=item.generate_cas_outputs_gcs_config_path,
        )
        _ = generate_cas_outputs_op()


#
#
# @dsl.pipeline(name="calculate_metrics", description="Calculate metrics for cas outputs")
def calculate_metrics_pipeline(calculate_metrics_config_paths: t.List[str]):
    with dsl.ParallelFor(calculate_metrics_config_paths) as item:
        calculate_metrics_op = create_job(
            component_func=job_components.benchmarking.calculate_metrics,
            component_name=constants.CALCULATE_METRICS_COMPONENT_NAME,
            gcs_config_path=item.calculate_metrics_gcs_config_path,
        )
        _ = calculate_metrics_op()


#
#
# @dsl.pipeline(name="bq_ops_create_avro_files", description="Create avro files for BigQuery ingest")
# def bq_ops_create_avro_files(pipeline_config_paths: t.List[str]):
#     """
#
#     :rtype: dsl.Pipeline
#     """
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         create_avro_files_op = create_job(
#             component_func=job_components.bq_ops.create_avro_files,
#             component_name=constants.BQ_OPS_CREATE_AVRO_FILES_COMPONENT_NAME,
#             gcs_config_path=item.bq_ops_create_avro_files_gcs_config_path,
#         )
#         _ = create_avro_files_op()
#
#
# @dsl.pipeline(name="bq_ops_ingest_data", description="Ingest data to BigQuery")
# def bq_ops_ingest_data(pipeline_config_paths: t.List[str]):
#     """
#
#     :rtype: dsl.Pipeline
#     """
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         ingest_files_op = create_job(
#             component_func=job_components.bq_ops.ingest_data,
#             component_name=constants.BQ_OPS_INGEST_DATA_COMPONENT_NAME,
#             gcs_config_path=item.bq_ops_ingest_data_gcs_config_path,
#         )
#         _ = ingest_files_op()
#
#
# @dsl.pipeline(
#     name="bq_ops_create_avro_files_and_ingest_data", description="Create avro files and ingest them into BigQuery"
# )
# def bq_ops_create_avro_files_and_ingest_data(
#     bq_ops_create_avro_files_config_paths: t.List[str], bq_ops_ingest_data_config_paths: t.List[str]
# ):
#     bq_ops_create_avro_files_task = bq_ops_create_avro_files(
#         pipeline_config_paths=bq_ops_create_avro_files_config_paths
#     )
#     bq_ops_ingest_data_task = bq_ops_ingest_data(pipeline_config_paths=bq_ops_ingest_data_config_paths)
#     bq_ops_ingest_data_task.after(bq_ops_create_avro_files_task)
#
#
# @dsl.pipeline(name="bq_ops_precalculate_fields", description="Precalculate fields in BigQuery")
# def bq_ops_precalculate_fields(pipeline_config_paths: t.List[str]):
#     with dsl.ParallelFor(pipeline_config_paths) as item:
#         precalculate_fields_op = create_job(
#             component_func=job_components.bq_ops.precalculate_fields,
#             component_name=constants.BQ_OPS_PRECALCULATE_FIELDS_COMPONENT_NAME,
#             gcs_config_path=item.bq_ops_precalculate_fields_gcs_config_path,
#         )
#         _ = precalculate_fields_op()
#
#
@dsl.pipeline(name="bq_ops_prepare_extract", description="Prepare tables for extract in BigQuery")
def bq_ops_prepare_extract(pipeline_config_paths: t.List[str]):
    """

    :rtype: dsl.Pipeline
    """
    with dsl.ParallelFor(pipeline_config_paths) as item:
        prepare_extract_op = create_job(
            component_func=job_components.bq_ops.prepare_extract,
            component_name=constants.BQ_OPS_PREPARE_EXTRACT_COMPONENT_NAME,
            gcs_config_path=item.bq_ops_prepare_extract_gcs_config_path,
        )
        _ = prepare_extract_op()


@dsl.pipeline(name="bq_ops_extract", description="Extract data from BigQuery")
def bq_ops_extract(pipeline_config_paths: t.List[str]):
    """

    :rtype: dls.Pipeline
    """
    with dsl.ParallelFor(pipeline_config_paths) as item:
        extract_op = create_job(
            component_func=job_components.bq_ops.extract,
            component_name=constants.BQ_OPS_EXTRACT_COMPONENT_NAME,
            gcs_config_path=item.bq_ops_extract_gcs_config_path,
        )
        _ = extract_op()


@dsl.pipeline(
    name="bq_ops_prepare_and_extract", description="Prepare extract tables in BigQuery and extract data to GCS bucket"
)
def bq_ops_prepare_and_extract(prepare_extract_config_paths: t.List[str], extract_config_paths: t.List[str]) -> None:
    prepare_extract_task = bq_ops_prepare_extract(pipeline_config_paths=prepare_extract_config_paths)
    extract_task = bq_ops_extract(pipeline_config_paths=extract_config_paths)
    extract_task.after(prepare_extract_task)

    # with dsl.ParallelFor(prepare_extract_config_paths) as item:
    #     prepare_extract_op = create_job(
    #         component_func=job_components.bq_ops.prepare_extract,
    #         component_name=constants.BQ_OPS_PREPARE_EXTRACT_COMPONENT_NAME,
    #         gcs_config_path=item.bq_ops_prepare_extract_gcs_config_path,
    #     )
    # _ = prepare_extract_op()
    #
    # with dsl.ParallelFor(extract_config_paths) as item:
    #     extract_op = create_job(
    #         component_func=job_components.bq_ops.extract,
    #         component_name=constants.BQ_OPS_EXTRACT_COMPONENT_NAME,
    #         gcs_config_path=item.bq_ops_extract_gcs_config_path,
    #     )
    # _ = extract_op()
    #
    # dsl
