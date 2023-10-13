task precalculate_fields {
    String version = "precalculate_fields_v1.0.0"
    String docker_image
    String dataset
    String fields

    command {
        cd /app
        echo "START INGESTING"
        python casp/services/bq_ops/ingest_files_to_bq/main.py \
        --dataset ${dataset} \
        --fields ${fields}
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

workflow CASPrecalculateFields {
    meta {
        description: "Precalculate Fields"
    }
    call precalculate_fields
}