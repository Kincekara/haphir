version 1.1

task estimate_genome_size {
    input {
        File long_fq
        Int cpu = 1
    }

    command <<<
        set -euo pipefail

        # version 
        lrge --version > VERSION

        # find genome size
        lrge -t ~{cpu} ~{long_fq} -o gsize.txt
    >>>

    output {
        String lrge_version = read_string("VERSION")
        String genome_size = read_string("gsize.txt")
    }

    runtime {
        container: "staphb/lrge:0.2.1"
        cpu: cpu
        memory: "1 GiB"
        preemptible: 0
        maxRetries: 3
    }
}