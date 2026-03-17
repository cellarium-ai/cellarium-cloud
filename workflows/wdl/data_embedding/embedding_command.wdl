task embed {
    String version = "train_v1.0.0"
    String docker_image
    String bucket_name
    String data_storage_path
    String dm_storage_path
    String output_storage_path
    String running_script

    command {
        echo $(pwd)
        cd /app
        python ${running_script} \
        --bucket_name=${bucket_name} \
        --data_storage_path=${data_storage_path} \
        --dm_storage_path=${dm_storage_path} \
        --output_storage_path=${output_storage_path} \
    }

    runtime {
        docker: docker_image
        bootDiskSizeGb: 100
        memory: "26G"
        cpu: 4
        gpuCount: 1
        gpuType:  "nvidia-tesla-t4"
        zones: "us-east1-d us-east1-c"
        maxRetries: 0
        preemptible_tries: 0
    }
}

workflow CASPEmbedData {
    String bucket_name
    String data_storage_path
    String dm_storage_path
    String output_storage_path
    String running_script

    meta {
        description: "Embedding Workflow"
    }

    call embed {
        input:
            bucket_name = bucket_name,
            data_storage_path = data_storage_path,
            dm_storage_path = dm_storage_path,
            output_storage_path = output_storage_path,
            running_script = running_script,
    }
}