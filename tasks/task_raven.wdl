version 1.0

task raven_asm {
    input {
        String id
        File long_fq
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version
        raven --version > VERSION

        # assemble with raven
        raven \
        --threads ~{cpu} \
        --kmer-len 29 \
        --window-len 9 \
        --graphical-fragment-assembly ~{id}.raven.gfa \
        ~{long_fq} > ~{id}.raven.fasta
    >>>

    output {
        String raven_version = read_string("VERSION")
        File assembly_graph = "~{id}.raven.gfa"
        File assembly_fasta = "~{id}.raven.fasta"
    }

    runtime {
        docker: "quay.io/biocontainers/raven-assembler:1.8.3--h5ca1c30_3"
        cpu: cpu
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}