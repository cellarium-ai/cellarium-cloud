task prepare_extract {
    String version = "prepare_extract_v1.0.0"
    String bq_dataset
    String extract_table_prefix

    command {
        cd /app
        echo "START PREPARING EXTRACT"
        python casp/services/bq_ops/prepare_extract/main.py \
        --dataset ${bq_dataset} \
        --extract_table_prefix ${extract_table_prefix}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:fg-50m-2"
        bootDiskSizeGb: 10
        memory: "2G"
        cpu: 1
        zones: "us-central1-a"
    }
}

workflow CASPrepareExtractBQ {

    meta {
        description: "Prepare Data for extract in BigQuery extract tables"
    }

    call prepare_extract
}