task extract {
    String version = "extract_v1.1.0"
    String project_id
    String bq_dataset
    String extract_table_prefix
    Int start_bin
    Int end_bin
    String output_bucket_name
    String output_bucket_directory

    command {
        cd /app
        echo "START EXTRACTING"
        python casp/services/bq_ops/extract/main.py \
        --dataset ${bq_dataset} \
        --extract_table_prefix ${extract_table_prefix} \
        --start_bin ${start_bin} \
        --end_bin ${end_bin} \
        --output_bucket_name ${output_bucket_name} \
        --output_bucket_directory ${output_bucket_directory}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:fg-50m-35"
        bootDiskSizeGb: 10
        memory: "20G"
        cpu: 10
        zones: "us-central1-a"
        preemptible: 3
    }
}

workflow CASExtractBQ {
    Array[Array[Int]] bin_borders

    meta {
        description: "Extract Data from BigQuery to anndata files"
    }

    scatter(bin_border in bin_borders) {
        call extract {
            input:
                start_bin = bin_border[0],
                end_bin = bin_border[1]
        }
    }
}