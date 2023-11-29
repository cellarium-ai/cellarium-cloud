task prepare_extract {
    String version = "prepare_extract_v1.0.0"
    String docker_image
    String bq_dataset
    String extract_table_prefix
    String fq_allowed_original_feature_ids
    Int extract_bin_size
    String bucket_name
    String extract_bucket_path
    String filter_by_organism
    Boolean filter_by_is_primary_data
    String? filter_by_diseases
    String? obs_columns_to_include

    command {
        cd /app
        echo "START PREPARING EXTRACT"
        python casp/services/bq_ops/prepare_extract/main.py \
        --dataset ${bq_dataset} \
        --extract_table_prefix ${extract_table_prefix} \
        --fq_allowed_original_feature_ids ${fq_allowed_original_feature_ids} \
        --extract_bin_size ${extract_bin_size} \
        --bucket_name ${bucket_name} \
        --extract_bucket_path ${extract_bucket_path} \
        --filter_by_organism "${filter_by_organism}" \
        --filter_by_is_primary_data ${filter_by_is_primary_data} \
        --filter_by_diseases "${filter_by_diseases}" \
        --obs_columns_to_include "${obs_columns_to_include}"
    }

    runtime {
        docker: docker_image
        bootDiskSizeGb: 20
        memory: "10G"
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