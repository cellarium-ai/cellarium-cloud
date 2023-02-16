task load_dataset {
    String version = "ingest_to_bq_v1.0.0"
    String bq_dataset
    String gcs_input_bucket
    String gcs_file_path
    String gcs_stage_prefix
    Int cas_cell_index_start
    Int cas_feature_index_start


    command {
        cd /app
        echo "START INGESTING"
        python casp/services/bq_ops/load_dataset/main.py \
        --dataset ${bq_dataset} \
        --gcs_input_bucket ${gcs_input_bucket} \
        --gcs_file_path ${gcs_file_path} \
        --gcs_stage_prefix ${gcs_stage_prefix} \
        --cas_cell_index_start ${cas_cell_index_start} \
        --cas_feature_index_start ${cas_feature_index_start}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-manual/cas-pytorch:fg-ingest-01"
        bootDiskSizeGb: 50
        memory: "4G"
        cpu: 1
        zones: "us-central1-a"
        maxRetries: 0
        preemptible_tries: 3
    }
}

workflow CASIngestData {
    Array[Object] convert_args

    meta {
        description: "Ingest data to BigQuery"
    }

    scatter(arg in convert_args) {
        call load_dataset {
            input:
                gcs_file_path = arg.df_filename,
                cas_cell_index_start = arg.cas_cell_index,
                cas_feature_index_start = arg.cas_feature_index
        }
    }
}