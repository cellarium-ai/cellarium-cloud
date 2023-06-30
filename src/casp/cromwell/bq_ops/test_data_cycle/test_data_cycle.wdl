task test_data_cycle {
    String version = "test_data_cycle.1.0"

    command {
        cd /app
        echo "START DATA CYCLE TEST"
        python casp/services/bq_ops/test_data_cycle/main.py
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/dsp-cell-annotation-service/cas-services-cicd/cas-pytorch:june-release-11"
        bootDiskSizeGb: 40
        memory: "30G"
        cpu: 18
        zones: "us-central1-a"
        preemptible: 3
    }
}

workflow CASTestDataCycle {
    meta {
        description: "Test CAS data pipeline"
    }
    call test_data_cycle
}