task ingest_files_to_bq {
    String version = "ingest_files_to_bq_v1.0.0"
    String docker_image
    String dataset
    String gcs_bucket_name
    String gcs_stage_dir

    command {
        cd /app
        echo "START INGESTING"
        python casp/services/bq_ops/ingest_files_to_bq/main.py \
        --dataset ${dataset} \
        --gcs_bucket_name ${gcs_bucket_name} \
        --gcs_stage_dir ${gcs_stage_dir}
    }

    runtime {
        docker: docker_image
        bootDiskSizeGb: 25
        memory: "1G"
        cpu: 1
        zones: "us-central1-a"
        maxRetries: 3
        preemptible_tries: 3
    }
}

workflow CASIngestFilesToBQ {
    meta {
        description: "Ingest files to BigQuery"
    }
    call ingest_files_to_bq
}