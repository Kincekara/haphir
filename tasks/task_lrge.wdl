version 1.0

task estimate_genome_size {
    input {
        File long_fq
        Int cpu = 4
    }

    command <<<
        set -euo pipefail

        # version 
        lrge --version > VERSION

        # find genome size
        lrge \
        -P pb \
        -t ~{cpu} \
        -o gsize.txt \
        ~{long_fq}   

        # round genome size
        printf "%dm\n" $(( ($(cat gsize.txt) +500000)/1000000 )) > GSIZE
    >>>

    output {
        String lrge_version = read_string("VERSION")
        File lrge_gs = "gsize.txt"
        String genome_size = read_string("GSIZE")
    }

    runtime {
        docker: "staphb/lrge:0.2.1"
        cpu: cpu
        memory: "8 GiB"
        disks: "local-disk 50 HDD"
        preemptible: 0
        maxRetries: 3
    }
}