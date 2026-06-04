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
        --identity 0.99 \
        --kmer-len 29 \
        --window-len 9 \
        --polishing-rounds 1 \
        --graphical-fragment-assembly ~{id}.raven.gfa \
        ~{long_fq} > ~{id}.raven.fasta

        # get contig lengths
        echo "Raven" > ~{id}.raven.ctg_len.txt
        awk -F'LN:i:' '/^>/{split($2,a," "); print a[1]}' ~{id}.raven.fasta | sort -nr >> ~{id}.raven.ctg_len.txt
    >>>

    output {
        String raven_version = read_string("VERSION")
        File assembly_graph = "~{id}.raven.gfa"
        File assembly_fasta = "~{id}.raven.fasta"
        File ctg_len = "~{id}.raven.ctg_len.txt"
    }

    runtime {
        docker: "staphb/raven:1.8.3-noble"
        cpu: cpu
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}