task train {
    String version = "train_v1.0.0"
    String training_script
    String bucket_name
    String data_storage_path
    String checkpoint_save_path
    Int n_components
    Int q
    Int batch_size
    Boolean use_gpu

    command {
        echo $(pwd)
        cd /app
        echo $training_script
        python casp/ml/services/training/pca_lowrank/train.py \
        --bucket_name=${bucket_name} \
        --data_storage_path=${data_storage_path} \
        --checkpoint_save_path=${checkpoint_save_path} \
        --n_components=${n_components} \
        --q=${q} \
        --batch_size=${batch_size} \
        --use_gpu=${use_gpu}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch-cuda:fg-pca-lowrank-01"
        bootDiskSizeGb: 50
        memory: "26G"
        cpu: 4
        gpuCount: 1
        gpuType:  "nvidia-tesla-t4"
        zones: "us-east1-d us-east1-c"
        maxRetries: 0
        preemptible_tries: 0
    }
}

workflow CASTrainPCALowrank {
    String training_script
    String bucket_name
    String data_storage_path
    String checkpoint_save_path
    Int n_components
    Int q
    Int batch_size
    Boolean use_gpu

    meta {
        description: "Train workflow"
    }

    call train {
        input:
            training_script = training_script,
            bucket_name = bucket_name,
            data_storage_path = data_storage_path,
            checkpoint_save_path = checkpoint_save_path,
            n_components = n_components,
            q = q,
            batch_size = batch_size,
            use_gpu = use_gpu
    }
}