task precalculate_fields {
    String version = "precalculate_fields_v1.0.0"
    String docker_image
    String dataset
    String fields

    command {
        cd /app
        echo "START INGESTING"
        python casp/services/bq_ops/precalculate_fields/main.py \
        --dataset ${dataset} \
        --fields_str ${fields}
    }

    runtime {
        docker: docker_image
        bootDiskSizeGb: 25
        memory: "4G"
        cpu: 1
        zones: "us-central1-a"
        maxRetries: 3
    }
}

workflow CASPrecalculateFields {
    meta {
        description: "Precalculate Fields"
    }
    call precalculate_fields
}