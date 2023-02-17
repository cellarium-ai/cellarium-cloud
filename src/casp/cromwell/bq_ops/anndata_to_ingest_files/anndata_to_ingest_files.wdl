task anndata_to_ingest_files {
    String version = "ingest_to_bq_v1.0.0"
    String bq_dataset
    String gcs_input_bucket
    String gcs_file_path
    String gcs_ingest_path
    Int cas_cell_index_start
    Int cas_feature_index_start


    command {
        cd /app
        echo "START INGESTING"
        python casp/services/bq_ops/anndata_to_ingest_files/main.py \
        --dataset ${bq_dataset} \
        --gcs_input_bucket ${gcs_input_bucket} \
        --gcs_file_path ${gcs_file_path} \
        --gcs_ingest_path ${gcs_ingest_path} \
        --cas_cell_index_start ${cas_cell_index_start} \
        --cas_feature_index_start ${cas_feature_index_start}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-manual/cas-pytorch:fg-ingest-01"
        bootDiskSizeGb: 50
        memory: "16G"
        cpu: 4
        zones: "us-central1-a"
        maxRetries: 3
        preemptible_tries: 3
    }
}

workflow CASAnndataToIngestFiles {
    Array[Object] convert_args

    meta {
        description: "Create avro "
    }

    scatter(arg in convert_args) {
        call anndata_to_ingest_files {
            input:
                gcs_file_path = arg.df_filename,
                cas_cell_index_start = arg.cas_cell_index,
                cas_feature_index_start = arg.cas_feature_index
        }
    }
}