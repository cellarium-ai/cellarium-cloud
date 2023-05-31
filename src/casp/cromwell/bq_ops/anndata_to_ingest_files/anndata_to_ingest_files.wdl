task anndata_to_ingest_files {
    String version = "anndata_to_files_v1.0.0"
    String gcs_input_bucket
    String gcs_file_path
    String gcs_stage_dir
    Int cas_cell_index_start
    Int cas_feature_index_start
    String original_feature_id_lookup
    Boolean load_uns_data
    String uns_meta_keys


    command {
        cd /app
        echo "START MAKING INGEST FILES"
        python casp/services/bq_ops/anndata_to_ingest_files/main.py \
        --gcs_input_bucket ${gcs_input_bucket} \
        --gcs_file_path ${gcs_file_path} \
        --gcs_stage_dir ${gcs_stage_dir} \
        --cas_cell_index_start ${cas_cell_index_start} \
        --cas_feature_index_start ${cas_feature_index_start} \
        --original_feature_id_lookup ${original_feature_id_lookup} \
        --load_uns_data ${load_uns_data} \
        --uns_meta_keys ${uns_meta_keys}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:fg-50m-33"
        bootDiskSizeGb: 50
        memory: "64G"
        cpu: 32
        zones: "us-central1-a"
        maxRetries: 3
        preemptible_tries: 3
    }
}

workflow CASAnndataToIngestFiles {
    Array[Object] convert_args

    meta {
        description: "Create ingest files"
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