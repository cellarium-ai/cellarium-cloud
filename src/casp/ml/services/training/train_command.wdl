task train {
    String version = "train_v1.0.0"
    String bucket_name
    String data_storage_path
    String training_script

    command {
        echo $(pwd)
        cd /app
        python ${training_script} \
        --bucket_name=${bucket_name} \
        --data_storage_path=${data_storage_path}
    }

    runtime {
        docker: "us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_training-gpu:1.0"
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

workflow CASPModelTrain {
    String bucket_name
    String data_storage_path
    String training_script

    meta {
        description: "Train workflow"
    }

    call train {
        input:
            bucket_name = bucket_name,
            data_storage_path = data_storage_path,
            training_script = training_script,
    }
}