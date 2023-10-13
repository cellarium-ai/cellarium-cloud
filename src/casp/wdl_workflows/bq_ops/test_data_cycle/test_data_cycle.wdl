task test_data_cycle {
    String version = "test_data_cycle.1.0"
    String docker_image

    command {
        cd /app
        echo "START DATA CYCLE TEST"
        /home/user/micromamba/bin/pytest tests/infrastructure --log-cli-level=INFO
    }

    runtime {
        docker: docker_image
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