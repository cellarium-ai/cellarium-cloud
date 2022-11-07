task train {

    String hello_input
    String version = "train_v1.0.0"
    String python_environment

    command {
        echo ${hello_input}
        echo $(pwd)
        cd /app
        python /app/casp/ml/training/pca/train.py
        python -u<<CODE
        chars = """Hello World!"""
        print(chars)
        import sys
        print(sys.version)
        CODE
    }

    runtime {
        docker: "us-east4-docker.pkg.dev/dsp-cell-annotation-service/casp-pca/casp_pca_training:1.0"
        bootDiskSizeGb: 100
        memory: "26G"
        cpu: 4
        zones: "us-east1-d us-east1-c"
        maxRetries: 0
        preemptible_tries: 0
    }
}

workflow CASPModelTrain {
    String python_environment

    meta {
        description: "Train workflow"
    }

    call train {
        input:
            hello_input = "Hello World",
            python_environment = python_environment
    }
}