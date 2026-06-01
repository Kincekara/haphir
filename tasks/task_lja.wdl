version 1.0

task lja_asm {
    input {
        String id
        File long_fq
        Int cpu = 8
    }

    command <<<
        set -euo pipefail
    
        # assemble with lja
        lja \
        -t ~{cpu} \
        -o out \
        --reads ~{long_fq} > lja.out.txt 2> lja.err.txt

        # rename output
        mv out/assembly.fasta ~{id}.lja.fasta
        mv out/mdbg.gfa ~{id}.lja.gfa

    >>>

    output {
        String lja_version = "0.2"
        File assembly_fasta = "~{id}.lja.fasta"
        File assembly_graph = "~{id}.lja.gfa"
    }

    runtime {
        docker: "kincekara/lja:0.2"
        cpu: cpu
        memory: "32 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}