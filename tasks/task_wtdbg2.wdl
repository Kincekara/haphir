version 1.0

task wtdbg2_asm {
    input {
        String id
        File long_fq
        String genome_size
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        wtdbg2 --version > VERSION

        # assemble with wtdb2
        wtdbg2 \
        -x ccs \
        -t ~{cpu} \
        -i ~{long_fq} \
        -g ~{genome_size} \
        -o ~{id}

        # derive consensus
        wtpoa-cns \
        -t ~{cpu} \
        -i ~{id}.ctg.lay.gz -fo ~{id}.wtdbg2.fasta
    >>>

    output {
        String wtdbg2_version = read_string("VERSION")
        File assembly_fasta = "~{id}.wtdbg2.fasta"
    }

    runtime {
        docker: "staphb/wtdbg2:2.5"
        cpu: cpu
        memory: "8 GiB"
        preemptible: 2
        maxRetries: 5
    }
}

